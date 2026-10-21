#!/usr/bin/env python3
"""VCF-style variant table (no header: chr pos ref alt [extra columns...]) -> VCF.

Usage: vcftab2vcf.py <table> <ref.fa.gz> <out.vcf>
ID = row number (r<N>); INFO/ORIG = original chr:pos:ref:alt; INFO/X = extra columns joined with '|';
INFO/REFOK = 1 if REF matches the genome. Rows on contigs missing from the genome are skipped (counted).
"""
import sys, pysam

tab, fa, out = sys.argv[1:4]
g = pysam.FastaFile(fa)
contigs = set(g.references)
st = dict(rows=0, refmismatch=0, nocontig=0)
with open(tab) as fh, open(out, "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    for c in g.references:
        if "_" not in c:
            o.write(f"##contig=<ID={c},length={g.get_reference_length(c)}>\n")
    o.write('##INFO=<ID=ORIG,Number=1,Type=String,Description="original chr:pos:ref:alt">\n'
            '##INFO=<ID=X,Number=1,Type=String,Description="extra columns, |-joined">\n'
            '##INFO=<ID=REFOK,Number=1,Type=Integer,Description="REF matches genome">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for n, l in enumerate(fh, 1):
        f = l.rstrip("\r\n").split("\t")
        if len(f) < 4 or not f[1].isdigit():
            continue
        st["rows"] += 1
        c, p, r, a = f[0], int(f[1]), f[2].upper(), f[3].upper()
        if c not in contigs:
            st["nocontig"] += 1
            continue
        ok = int(g.fetch(c, p - 1, p - 1 + len(r)).upper() == r)
        st["refmismatch"] += 1 - ok
        o.write(f"{c}\t{p}\tr{n}\t{r}\t{a}\t.\t.\tORIG={f[0]}:{f[1]}:{f[2]}:{f[3]};X={'|'.join(f[4:])};REFOK={ok}\n")
print(tab, st, file=sys.stderr)
