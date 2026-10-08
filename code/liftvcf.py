#!/usr/bin/env python3
"""Lift a (normalized) VCF to another build with the UCSC liftOver binary + chain.

Usage: liftvcf.py <in.vcf> <chain.gz> <target.fa.gz> <out.vcf> <out.unmap>

For every record the REF span [POS, POS+len(REF)) is lifted with `liftOver -minMatch=1`.
 - unmapped by liftOver               -> unmap, reason from liftOver (Deleted/Partially/Split/Duplicated)
 - lands on a non-primary contig       -> unmap "alt_contig"
 - lifted span length != len(REF)      -> unmap "span_length_changed"
 - minus strand: REF and ALT are reverse-complemented (re-anchoring is left to bcftools norm)
 - the .unmap file has LIFTED_AT (chr:start:end:strand[:target REF]) when liftOver did map the span
 - REF is always re-read from the target genome; if it differs from the original
   (strand-adjusted) REF the record gets REFDIFF=<old>; if the new REF equals the ALT,
   the record is unmapped as "ref_alt_swap" (genotype/AF would have to be inverted).
INFO STRAND=+/- is added.
"""
import sys, subprocess, tempfile, os, pysam

vcf, chain, fa, out, unmap = sys.argv[1:6]
genome = pysam.FastaFile(fa)
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(comp)[::-1]

header, recs = [], {}
with open(vcf) as fh:
    for l in fh:
        if l[0] == "#":
            if not l.startswith("##contig"):
                header.append(l)
            continue
        f = l.rstrip("\n").split("\t")
        recs[f[2]] = f

tmp = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(out)))
bed, lo, un = (os.path.join(tmp, x) for x in ("in.bed", "out.bed", "un.bed"))
with open(bed, "w") as b:
    for vid, f in recs.items():
        p = int(f[1]) - 1
        b.write(f"{f[0]}\t{p}\t{p+len(f[3])}\t{vid}\t0\t+\n")
subprocess.run(["liftOver", "-minMatch=1", bed, chain, lo, un], check=True, stderr=subprocess.DEVNULL)

reason, lifted_at = {}, {}
with open(un) as fh:
    why = "unmapped"
    for l in fh:
        if l[0] == "#":
            why = l[1:].strip().replace(" ", "_")
        else:
            reason[l.split("\t")[3]] = why

nstat = dict(mapped=0, refdiff=0, minus=0)
with open(out, "w") as o, open(unmap, "w") as u:
    hdr = "".join(header[:-1])
    o.write(hdr)
    for c in genome.references:
        if "_" not in c:
            o.write(f"##contig=<ID={c},length={genome.get_reference_length(c)}>\n")
    o.write('##INFO=<ID=STRAND,Number=1,Type=String,Description="chain strand">\n')
    o.write('##INFO=<ID=REFDIFF,Number=1,Type=String,Description="REF before liftover differed from target genome">\n')
    o.write(header[-1])
    u.write(header[-1].rstrip("\n") + "\tLIFTED_AT\tREASON\n")
    with open(lo) as fh:
        for l in fh:
            c, s, e, vid, _, strand = l.rstrip("\n").split("\t")[:6]
            f = recs[vid]
            s, e = int(s), int(e)
            ref, alt = f[3].upper(), f[4].upper()
            lifted_at[vid] = f"{c}:{s+1}:{e}:{strand}"
            if "_" in c or c not in genome.references:
                reason[vid] = "alt_contig:" + c; continue
            if e - s != len(ref):
                reason[vid] = "span_length_changed"; continue
            if strand == "-":
                ref, alt = rc(ref), rc(alt)
                nstat["minus"] += 1
            newref = genome.fetch(c, s, e).upper()
            lifted_at[vid] += f":{newref}"
            info = f[7] + f";STRAND={strand}"
            if newref != ref:
                if newref == alt:
                    reason[vid] = "ref_alt_swap"; continue
                info += f";REFDIFF={ref}"
                nstat["refdiff"] += 1
            nstat["mapped"] += 1
            o.write("\t".join([c, str(s + 1), vid, newref, alt] + f[5:7] + [info] + f[8:]) + "\n")
    for vid, why in reason.items():
        u.write("\t".join(recs[vid]) + "\t" + lifted_at.get(vid, ".") + "\t" + why + "\n")
print(os.path.basename(vcf), "->", os.path.basename(out), nstat, "unmapped:", len(reason), file=sys.stderr)
subprocess.run(["rm", "-r", tmp])
