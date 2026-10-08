#!/usr/bin/env python3
"""Convert an Iranome ANNOVAR-style table (Chr Start End Ref Alt AF) to VCF.

Usage: annovar2vcf.py <table> <ref.fa(.gz, faidx'd)> <out.vcf> <tag> [--anchored-ins]

- SNV / MNP / complex: POS=Start, REF=Ref, ALT=Alt
- deletion (Alt '-'):  POS=Start-1, REF=anchor+Ref, ALT=anchor
- insertion (Ref '-'): POS=Start,   REF=base@Start, ALT=base+Alt   (ANNOVAR convention)
Every record gets ID=<tag>_<line number> so it can be traced through liftover.
INFO: AF, ORIG=chr:start:ref:alt (original row), REFOK=1/0 (table Ref matches genome).
"""
import sys, pysam

tab, fa, out, tag = sys.argv[1:5]
genome = pysam.FastaFile(fa)
contigs = set(genome.references)
stats = dict(total=0, refmismatch=0, nocontig=0)

def chrom_name(c):
    c = c if c.startswith("chr") else "chr" + c
    return "chrM" if c == "chrMT" else c

extra = []
with open(tab) as fh, open(out, "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    for c in genome.references:
        if "_" not in c:
            o.write(f"##contig=<ID={c},length={genome.get_reference_length(c)}>\n")
    o.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Iranome AF">\n')
    o.write('##INFO=<ID=ORIG,Number=1,Type=String,Description="Original table row chr:start:end:ref:alt">\n')
    o.write('##INFO=<ID=REFOK,Number=1,Type=Integer,Description="Table Ref matches genome">\n')
    first = fh.readline()
    extra = first.rstrip("\n").lstrip("#").split("\t")[6:] if first.startswith("#") else []
    for x in extra:
        o.write(f'##INFO=<ID={x},Number=1,Type=Integer,Description="{x} from input table">\n')
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    fh.seek(0)
    for n, line in enumerate(fh, 1):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        c, s, e, r, a, af = cols[:6]
        xinfo = "".join(f";{k}={v}" for k, v in zip(extra, cols[6:]))
        stats["total"] += 1
        ch = chrom_name(c)
        if ch not in contigs:
            stats["nocontig"] += 1
            continue
        s, e = int(s), int(e)
        if a == "-":                         # deletion
            pos = s - 1
            anc = genome.fetch(ch, pos - 1, pos).upper()
            g = genome.fetch(ch, s - 1, s - 1 + len(r)).upper()
            ref, alt = anc + r, anc
        elif r == "-":                       # insertion after Start
            pos = s
            anc = genome.fetch(ch, pos - 1, pos).upper()
            g = ""
            ref, alt = anc, anc + a
            r = ""
        else:
            pos = s
            g = genome.fetch(ch, s - 1, s - 1 + len(r)).upper()
            ref, alt = r, a
        ok = int(g == r.upper())
        if not ok:
            stats["refmismatch"] += 1
        o.write(f"{ch}\t{pos}\t{tag}_{n}\t{ref}\t{alt}\t.\t.\t"
                f"AF={af}{xinfo};ORIG={c}:{s}:{e}:{r or '-'}:{a};REFOK={ok}\n")
print(tab, stats, file=sys.stderr)
