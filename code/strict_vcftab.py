#!/usr/bin/env python3
"""Strict hg38 set for a VCF-style table (chr pos ref alt [extra...], no header): keep a variant only if
every check agrees (same checks as strict_filter.py), otherwise omit it.

Usage: strict_vcftab.py <table> <work dir> <transanno dir> <out prefix> [--flank 25]
  work dir:      output of vcftab_roundtrip.sh (orig19.norm.vcf, fwd38.norm.vcf, rt.roundtrip_all.tsv)
  transanno dir: contains ucsc_chain/ and minimap2_chain/ from transanno_roundtrip.sh
Checks: UCSC round trip identical; transanno (UCSC chain) and transanno (minimap2 chain) round trip identical,
no REF/ALT swap and the same hg38 position/alleles; no hg38 reference change; <flank> bp either side identical.
Output
  <prefix>_hg38_strict.txt        kept variants, same structure as the input, hg38 (left-aligned), sorted
  <prefix>_not_lifted_hg19.txt    omitted variants: the original input rows, unchanged, in input order
  <prefix>_not_lifted_reasons.tsv omitted variants with the first failed check and the UCSC hg38 position
  <prefix>_strict_stats.tsv       per chromosome: total, kept, omitted, % omitted; then counts per reason
"""
import sys, collections, pysam

table, work, ta, prefix = sys.argv[1:5]
F = int(sys.argv[sys.argv.index("--flank") + 1]) if "--flank" in sys.argv else 25
R = "/root/naser/liftover/ref"
g19, g38 = pysam.FastaFile(f"{R}/hg19.fa"), pysam.FastaFile(f"{R}/hg38.fa")
cache = {}
def seq(g, tag, c):
    if (tag, c) not in cache:
        cache[(tag, c)] = g.fetch(c).upper()
    return cache[(tag, c)]
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def rt(p, swaps=None):
    sw = set(open(swaps).read().split()) if swaps else set()
    d = {}
    with open(p) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 7)
            d[f[0]] = ("SWAPPED" if f[0] in sw else f[1], f[5])
    return d


U = rt(f"{work}/rt.roundtrip_all.tsv")
T = rt(f"{ta}/ucsc_chain/rt.roundtrip_all.tsv", f"{ta}/ucsc_chain/fwd_ref_changed.ids")
M = rt(f"{ta}/minimap2_chain/rt.roundtrip_all.tsv", f"{ta}/minimap2_chain/fwd_ref_changed.ids")
orig19, fwd = {}, {}
for l in open(f"{work}/orig19.norm.vcf"):
    if l[0] != "#":
        f = l.split("\t", 5)
        orig19[f[2]] = (f[0], int(f[1]), f[3])
for l in open(f"{work}/fwd38.norm.vcf"):
    if l[0] != "#":
        f = l.rstrip("\n").split("\t")
        fwd[f[2]] = (f[0], int(f[1]), f[3], f[4], "REFDIFF=" in f[7], "STRAND=-" in f[7])


def verdict(vid):
    us, uh = U.get(vid, ("missing", ""))
    if us != "IDENTICAL":
        return "UCSC round trip not identical", uh
    for name, D in (("transanno(UCSC chain)", T), ("transanno(minimap2 chain)", M)):
        s, h = D.get(vid, ("missing", ""))
        if s != "IDENTICAL" or h != uh:
            return f"{name}: " + ("REF/ALT swapped" if s == "SWAPPED" else "different hg38" if s == "IDENTICAL" else s), uh
    c38, p38, r38, a38, refdiff, minus = fwd[vid]
    if refdiff:
        return "hg38 reference allele differs from hg19", uh
    c, p, r = orig19[vid]
    s19, s38 = seq(g19, 19, c), seq(g38, 38, c38)
    if p - 1 - F < 0 or p38 - 1 - F < 0:
        return f"within {F} bp of a chromosome end", uh
    w19 = s19[p - 1 - F:p - 1 + len(r) + F]
    w38 = s38[p38 - 1 - F:p38 - 1 + len(r38) + F]
    if w38 != (rc(w19) if minus else w19):
        return f"flanking {F} bp differ between hg19 and hg38", uh
    return "", uh


ck = lambda c: (0, int(c[3:]), c) if c[3:].isdigit() else (1, {"X": 23, "Y": 24, "M": 25}.get(c[3:], 99), c)
stats = collections.defaultdict(collections.Counter)
reasons = collections.Counter()
kept = []
ids = sorted(U, key=lambda v: (ck(orig19[v][0]) if v in orig19 else (2, 0, ""), int(v[1:])))  # by chromosome for the cache
with open(table) as fh:
    rows = {f"r{n}": l.rstrip("\r\n") for n, l in enumerate(fh, 1)}
omitted = {}
last = None
for vid in ids:
    c = orig19[vid][0]
    if c != last:                      # keep only the current chromosome in memory
        cache.clear(); last = c
    why, uh = verdict(vid)
    stats[c]["total"] += 1
    if why:
        omitted[vid] = (why, uh)
        stats[c]["omitted"] += 1
        reasons[why] += 1
    else:
        x = rows[vid].split("\t")
        c38, p38, r38, a38 = fwd[vid][:4]
        kept.append((ck(c38), p38, r38, a38, "\t".join([c38, str(p38), r38, a38] + x[4:])))
        stats[c]["kept"] += 1

kept.sort(key=lambda k: k[:4])
with open(f"{prefix}_hg38_strict.txt", "w") as o:
    for k in kept:
        o.write(k[4] + "\n")
order = sorted(omitted, key=lambda v: int(v[1:]))
with open(f"{prefix}_not_lifted_hg19.txt", "w") as o, open(f"{prefix}_not_lifted_reasons.tsv", "w") as r:
    ncol = len(rows[order[0]].split("\t")) if order else 4
    r.write("\t".join(["chr", "pos", "ref", "alt"] + [f"col{i}" for i in range(5, ncol + 1)] + ["omitted_because", "UCSC_hg38"]) + "\n")
    for v in order:
        o.write(rows[v] + "\n")
        r.write(rows[v] + f"\t{omitted[v][0]}\t{omitted[v][1] or '.'}\n")
with open(f"{prefix}_strict_stats.tsv", "w") as s:
    s.write("chromosome\ttotal\tkept\tomitted\tpct_omitted\n")
    tot = collections.Counter()
    for c in sorted(stats, key=ck):
        v = stats[c]; tot.update(v)
        s.write(f"{c}\t{v['total']}\t{v['kept']}\t{v['omitted']}\t{100*v['omitted']/v['total']:.3f}\n")
    s.write(f"ALL\t{tot['total']}\t{tot['kept']}\t{tot['omitted']}\t{100*tot['omitted']/tot['total']:.3f}\n\n#reason\tcount\n")
    for w, v in reasons.most_common():
        s.write(f"#{w}\t{v}\n")
print(open(f"{prefix}_strict_stats.tsv").read())
