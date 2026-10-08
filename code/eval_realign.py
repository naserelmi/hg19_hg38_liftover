#!/usr/bin/env python3
"""Compare where reads land after realignment to hg38 with the chain-lifted hg38 position.

Usage: eval_realign.py <unmatched_variants_hg19.tsv> <sample hg19 BAM/CRAM> <sample .hg38.bam> <hg19.fa> <out.tsv>

For every variant, the primary reads covering its hg19 position are looked up (by name and mate) in the
hg38 BAM. Reads with MAPQ >= 20 are confident placements.
  status  no_reads           no read covers the variant in this sample
          ambiguous          < 50% of the reads are placed with MAPQ >= 20
          supports_liftover  >= 70% of confident reads cover the chain-lifted hg38 position
          elsewhere          most confident reads cover another hg38 locus (given in hg38_by_reads)
          mixed              confident reads split between loci
          no_liftover_pos    variant has no hg38 position (deleted); hg38_by_reads says where the reads went
"""
import sys, collections, pysam

vfile, bam19, bam38, ref19, out = sys.argv[1:6]
MINQ = 20

# hg38 placements of every realigned read: (name, is_read1) -> (chrom, start, end, mapq)
pl = {}
with pysam.AlignmentFile(bam38) as b:
    for r in b.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:
            continue
        pl[(r.query_name, r.is_read1)] = None if r.is_unmapped else (r.reference_name, r.reference_start, r.reference_end, r.mapping_quality)

cnt = collections.Counter()
bycat = collections.defaultdict(collections.Counter)
with pysam.AlignmentFile(bam19, reference_filename=ref19) as b, open(vfile) as fh, open(out, "w") as o:
    hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
    ix = {h: n for n, h in enumerate(hdr)}
    o.write("\t".join(hdr) + "\treads\tconfident\tat_lifted_pos\tstatus\thg38_by_reads\n")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        c, p = f[ix["hg19_Chr"]], int(f[ix["hg19_Pos"]])
        c38, p38 = f[ix["hg38_Chr"]], f[ix["hg38_Pos"]]
        reads = [r for r in b.fetch(c, p - 1, p)
                 if not (r.is_secondary or r.is_supplementary or r.is_unmapped)
                 and r.reference_start < p and r.reference_end and r.reference_end >= p]
        loc = [pl.get((r.query_name, r.is_read1)) for r in reads]
        conf = [x for x in loc if x and x[3] >= MINQ]
        at = sum(1 for x in conf if p38 != "." and x[0] == c38 and x[1] < int(p38) + 1 and x[2] >= int(p38) - 1)
        bins = collections.Counter(f"{x[0]}:{x[1] // 1000}kb" for x in conf)
        top = ",".join(f"{k}({v})" for k, v in bins.most_common(3))
        if not reads:
            st = "no_reads"
        elif len(conf) < 0.5 * len(reads):
            st = "ambiguous"
        elif p38 == ".":
            st = "no_liftover_pos"
        elif at >= 0.7 * len(conf):
            st = "supports_liftover"
        elif bins.most_common(1)[0][1] >= 0.7 * len(conf):
            st = "elsewhere"
        else:
            st = "mixed"
        cnt[st] += 1
        bycat[f[ix["category"]]][st] += 1
        o.write(l.rstrip("\n") + f"\t{len(reads)}\t{len(conf)}\t{at}\t{st}\t{top}\n")

sts = ["no_reads", "ambiguous", "supports_liftover", "elsewhere", "mixed", "no_liftover_pos"]
print("category\t" + "\t".join(sts))
for k, v in sorted(bycat.items()):
    print(k + "\t" + "\t".join(str(v[s]) for s in sts))
print("TOTAL\t" + "\t".join(str(cnt[s]) for s in sts))
