#!/usr/bin/env python3
"""Exon BED (dam_II.bed layout) from a refGene GTF.

Usage: gtf2exon_bed.py <refGene.gtf[.gz]> <out.bed> [--all-contigs] [--transcript-order]

Output columns: chrom, start (0-based), end, <gene>_<n>, strand, 0
 - only 'exon' features are used;
 - per gene (gene_id) and chromosome, the unique exon intervals of all its transcripts are numbered
   1..n in genomic order (left to right), as in dam_II.bed; --transcript-order numbers them 5'->3'
   instead (reversed on the minus strand);
 - strand comes from the GTF (per exon: the strand of the transcripts that use it);
 - primary chromosomes only (chr1-22, X, Y, M) unless --all-contigs.
Sorted by chromosome and start.
"""
import sys, gzip, re, collections

args = [x for x in sys.argv[1:] if not x.startswith("--")]
gtf, out = args
ALL, TXORDER = "--all-contigs" in sys.argv, "--transcript-order" in sys.argv
PRIMARY = re.compile(r"^chr([0-9]+|X|Y|M)$")
gid = re.compile(r'gene_id "([^"]+)"')

ex = collections.defaultdict(set)          # (gene, chrom) -> {(start0, end)}
strand = collections.defaultdict(set)      # (gene, chrom, start0, end) -> strands of the transcripts using it
with (gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf)) as fh:
    for l in fh:
        if l[0] == "#":
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or not (ALL or PRIMARY.match(f[0])):
            continue
        g = gid.search(f[8]).group(1)
        iv = (int(f[3]) - 1, int(f[4]))
        ex[(g, f[0])].add(iv)
        strand[(g, f[0]) + iv].add(f[6])

ck = lambda c: (0, int(c[3:]), c) if c[3:].isdigit() else (1, {"X": 23, "Y": 24, "M": 25}.get(c[3:], 99), c)
rows, mixed = [], 0
for (g, c), iv in ex.items():
    iv = sorted(iv)
    st = {x for a, b in iv for x in strand[(g, c, a, b)]}
    if TXORDER and st == {"-"}:
        iv = iv[::-1]
    for n, (a, b) in enumerate(iv, 1):
        s = strand[(g, c, a, b)]
        mixed += len(s) > 1
        rows.append((c, a, b, f"{g}_{n}", s.copy().pop() if len(s) == 1 else "."))
rows.sort(key=lambda r: (ck(r[0]), r[1], r[2], r[3]))
with open(out, "w") as o:
    for c, a, b, name, s in rows:
        o.write(f"{c}\t{a}\t{b}\t{name}\t{s}\t0\n")
print(f"{out}: {len(rows)} exons, {len(ex)} gene/chromosome groups"
      + (f", {mixed} exons used on both strands (strand '.')" if mixed else ""), file=sys.stderr)
