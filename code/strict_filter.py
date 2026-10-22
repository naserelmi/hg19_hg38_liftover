#!/usr/bin/env python3
"""Strict hg38 set: keep a variant only if every check agrees; otherwise omit it.

Usage: strict_filter.py <out dir> [--flank 25]
Run from /root/naser/liftover after run_pipeline.sh and both transanno runs. A variant is kept when
  1. UCSC liftOver round trip is identical                       (work_data/rt.roundtrip_all.tsv)
  2. transanno + UCSC chain: round trip identical, no REF/ALT swap, same hg38 position/alleles
  3. transanno + minimap2 chain: the same                          (transanno/<run>/rt.roundtrip_all.tsv)
  4. the hg38 reference allele equals the hg19 reference allele (no REFDIFF)
  5. the sequence <flank> bp either side of the variant is identical in hg19 and hg38
     (reverse-complemented for minus-strand chains)
Writes <out>/1_matched_variants_hg38_strict.tsv (same columns as results/1_matched_variants_hg38.tsv),
<out>/omitted_variants.tsv (first failed check per variant) and <out>/strict_summary.txt.
"""
import sys, collections, pysam

out = sys.argv[1]
F = int(sys.argv[sys.argv.index("--flank") + 1]) if "--flank" in sys.argv else 25
g19, g38 = pysam.FastaFile("ref/hg19.fa"), pysam.FastaFile("ref/hg38.fa")
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def rt(p, swaps=None):
    sw = set(open(swaps).read().split()) if swaps else set()
    d = {}
    with open(p) as fh:
        next(fh)
        for l in fh:
            f = l.rstrip("\n").split("\t")
            d[f[0]] = ("SWAPPED" if f[0] in sw else f[1], f[5])
    return d


U = rt("work_data/rt.roundtrip_all.tsv")
T = rt("transanno/ucsc_chain/rt.roundtrip_all.tsv", "transanno/ucsc_chain/fwd_ref_changed.ids")
M = rt("transanno/minimap2_chain/rt.roundtrip_all.tsv", "transanno/minimap2_chain/fwd_ref_changed.ids")

orig, fwd = {}, {}
for l in open("work_data/orig19.norm.vcf"):
    if l[0] != "#":
        f = l.rstrip("\n").split("\t")
        orig[f[2]] = (f[0], int(f[1]), f[3], f[4], dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv))
for l in open("work_data/fwd38.norm.vcf"):
    if l[0] != "#":
        f = l.rstrip("\n").split("\t")
        fwd[f[2]] = (f[0], int(f[1]), f[3], f[4], f[7])


def annovar(c, p, ref, alt):
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt:
        return c, p + 1, ref[1:], "_"
    if len(alt) > 1 and len(ref) == 1 and alt[0] == ref:
        return c, p, "_", alt[1:]
    return c, p, ref, alt


omit = collections.Counter()
kept = 0
with open(f"{out}/1_matched_variants_hg38_strict.tsv", "w") as k, open(f"{out}/omitted_variants.tsv", "w") as o:
    k.write("#Chr\tPos\tRef\tAlt\tHet\tHom\tAF\tAC\tAN\thg19_Chr\thg19_Pos\thg19_Ref\thg19_Alt\tStrand\thg38_ref_differs_from_hg19\n")
    o.write("#hg19_Chr\thg19_Pos\thg19_Ref\thg19_Alt\tHet\tHom\tAF\tAC\tAN\tomitted_because\tUCSC_hg38\n")
    for vid, (c, p, ref, alt, info) in orig.items():
        us, uh = U[vid]
        if us != "IDENTICAL":
            why = "UCSC round trip not identical"
        elif T[vid][0] != "IDENTICAL" or T[vid][1] != uh:
            why = "transanno(UCSC chain): " + ("REF/ALT swapped" if T[vid][0] == "SWAPPED" else
                  "different hg38" if T[vid][0] == "IDENTICAL" else T[vid][0])
        elif M[vid][0] != "IDENTICAL" or M[vid][1] != uh:
            why = "transanno(minimap2 chain): " + ("REF/ALT swapped" if M[vid][0] == "SWAPPED" else
                  "different hg38" if M[vid][0] == "IDENTICAL" else M[vid][0])
        elif "REFDIFF=" in fwd[vid][4]:
            why = "hg38 reference allele differs from hg19"
        else:
            c38, p38, r38, a38, inf38 = fwd[vid]
            w19 = g19.fetch(c, p - 1 - F, p - 1 + len(ref) + F).upper()
            w38 = g38.fetch(c38, p38 - 1 - F, p38 - 1 + len(r38) + F).upper()
            minus = "STRAND=-" in inf38
            why = "" if (w38 == (rc(w19) if minus else w19)) else f"flanking {F} bp differ between hg19 and hg38"
        g = [info["HET"], info["HOM"], info["AF"], info["AC"], info["AN"]]
        c0, s0, _, r0, a0 = info["ORIG"].split(":")
        h19 = [f"chr{c0}", s0, r0.replace("-", "_"), a0.replace("-", "_")]
        if why:
            omit[why.split(":")[0] if why.startswith("transanno") and ":" not in why[10:] else why] += 1
            o.write("\t".join(h19 + g + [why, uh or "."]) + "\n")
            continue
        c38, p38, r38, a38, inf38 = fwd[vid]
        k.write("\t".join(map(str, list(annovar(c38, p38, r38, a38)) + g + h19 +
                              ["-" if "STRAND=-" in inf38 else "+", "no"])) + "\n")
        kept += 1

tot = len(orig)
with open(f"{out}/strict_summary.txt", "w") as s:
    s.write(f"Iranome hg19 variants: {tot}\nkept in strict hg38 set: {kept} ({100*kept/tot:.3f}%)\n"
            f"omitted: {tot-kept} ({100*(tot-kept)/tot:.3f}%), first failed check:\n")
    for w, v in omit.most_common():
        s.write(f"  {v:>8d}  {w}\n")
print(open(f"{out}/strict_summary.txt").read())
