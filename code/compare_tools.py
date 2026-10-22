#!/usr/bin/env python3
"""Compare round-trip results of two liftover runs on the same input (per variant ID).

Usage: compare_tools.py <name A> <A rt.roundtrip_all.tsv> <name B> <B rt.roundtrip_all.tsv> <B fwd_ref_changed.ids|-> <out prefix>
Status per variant: matched (round trip identical), or the failure class.
A variant counts as 'swapped' for B when B swapped REF/ALT on the way to hg38 (transanno REF_CHANGED);
UCSC runs mark those as FAIL_hg19_to_hg38 (ref_alt_swap) themselves.
Writes <prefix>_crosstab.tsv, <prefix>_discordant.tsv and prints a summary.
"""
import sys, collections

na, fa, nb, fb, swb, out = sys.argv[1:7]
swapped = set(open(swb).read().split()) if swb != "-" else set()


def load(p, sw=frozenset()):
    d = {}
    with open(p) as fh:
        next(fh)
        for l in fh:
            f = l.rstrip("\n").split("\t")
            st = f[1]
            if f[0] in sw:
                st = "SWAPPED_" + st
            elif st == "FAIL_hg19_to_hg38" and f[8].startswith("ref_alt_swap"):
                st = "SWAPPED_not_lifted"
            d[f[0]] = (st, f[2], f[5], f[8])      # status, type, hg38, note
    return d


A, B = load(fa), load(fb, swapped)
simple = lambda s: "matched" if s == "IDENTICAL" else ("swapped" if s.startswith("SWAPPED") else
         "not_liftable_to_hg38" if s == "FAIL_hg19_to_hg38" else "not_liftable_back" if s == "FAIL_hg38_to_hg19" else "roundtrip_changed")
ct = collections.Counter()
same38 = diff38 = 0
with open(f"{out}_discordant.tsv", "w") as w:
    w.write(f"ID\ttype\t{na}_status\t{na}_hg38\t{na}_note\t{nb}_status\t{nb}_hg38\t{nb}_note\n")
    for vid, (sa, t, ha, xa) in A.items():
        sb, _, hb, xb = B.get(vid, ("missing", t, "", ""))
        ka, kb = simple(sa), simple(sb)
        ct[(ka, kb)] += 1
        if ka == kb == "matched":
            if ha == hb:
                same38 += 1
                continue
            diff38 += 1
        elif ka == kb:
            continue
        w.write(f"{vid}\t{t}\t{sa}\t{ha}\t{xa}\t{sb}\t{hb}\t{xb}\n")

cats = ["matched", "swapped", "not_liftable_to_hg38", "roundtrip_changed", "not_liftable_back"]
tot = len(A)
with open(f"{out}_crosstab.tsv", "w") as w:
    w.write(f"{na} \\ {nb}\t" + "\t".join(cats) + "\ttotal\n")
    for a in cats:
        w.write(a + "\t" + "\t".join(str(ct[(a, b)]) for b in cats) + f"\t{sum(ct[(a, b)] for b in cats)}\n")
    w.write("total\t" + "\t".join(str(sum(ct[(a, b)] for a in cats)) for b in cats) + f"\t{tot}\n")
print(open(f"{out}_crosstab.tsv").read())
for n, D_ in ((na, A), (nb, B)):
    m = sum(1 for v in D_.values() if simple(v[0]) == "matched")
    print(f"{n}: matched {m} ({100*m/tot:.3f}%)")
print(f"matched by both: {same38 + diff38}; same hg38 position/alleles: {same38}; different hg38: {diff38}")
