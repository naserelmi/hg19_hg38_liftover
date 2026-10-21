#!/usr/bin/env python3
"""Write results of the round-trip liftover of a VCF-style table.

Usage: vcftab_outputs.py <work dir> <original table> <out prefix>
Reads <work>/rt.roundtrip_all.tsv, fwd38.norm.vcf, fwd38.unmap, back19.unmap.
  <prefix>_hg38.txt        exact round-trip matches, same structure as the input (chr pos ref alt extra...),
                           hg38 coordinates/alleles (left-aligned on hg38), sorted by hg38 position, no header
  <prefix>_unmatched.tsv   every other row: original columns + category, reason, hg38 position, hg19 after round trip
  <prefix>_stats.tsv       per chromosome and overall: total, matched, mismatched, %, and counts per category
"""
import sys, collections

work, table, prefix = sys.argv[1:4]
orig = {}
with open(table) as fh:
    for n, l in enumerate(fh, 1):
        f = l.rstrip("\r\n").split("\t")
        if len(f) >= 4 and f[1].isdigit():
            orig[f"r{n}"] = f


def unmap(p):
    d = {}
    for l in open(p):
        if l[0] != "#":
            f = l.rstrip("\n").split("\t")
            d[f[2]] = (f[-2], f[-1])      # LIFTED_AT, REASON
    return d


fu, bu = unmap(f"{work}/fwd38.unmap"), unmap(f"{work}/back19.unmap")
h38 = {}
for l in open(f"{work}/fwd38.norm.vcf"):
    if l[0] != "#":
        f = l.split("\t", 5)
        h38[f[2]] = (f[0], int(f[1]), f[3], f[4])

CATS = {"ref_alt_swap": "ref_alt_swapped", "Deleted_in_new": "deleted_in_hg38", "Partially_deleted_in_new": "deleted_in_hg38",
        "alt_contig": "hg38_non_primary_contig", "span_length_changed": "broken", "Split_in_new": "broken"}
order = ["ref_alt_swapped", "deleted_in_hg38", "hg38_non_primary_contig", "broken", "not_liftable_back",
         "chromosome_changed", "duplicated_region_position_changed"]
stats = collections.defaultdict(collections.Counter)
reasons = collections.Counter()
matched = []
with open(f"{work}/rt.roundtrip_all.tsv") as fh, open(f"{prefix}_unmatched.tsv", "w") as u:
    next(fh)
    ncol = max(len(v) for v in orig.values())
    u.write("\t".join(["chr", "pos", "ref", "alt"] + [f"col{i}" for i in range(5, ncol + 1)] +
                      ["category", "reason", "hg38", "hg19_after_roundtrip"]) + "\n")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        vid, st = f[0], f[1]
        o = orig[vid]
        stats[o[0]]["total"] += 1
        if st == "IDENTICAL":
            c, p, r, a = h38[vid]
            matched.append((c, p, r, a, o[4:]))
            stats[o[0]]["matched"] += 1
            continue
        if st == "FAIL_hg19_to_hg38":
            at, why = fu[vid]
            cat = CATS.get(why.split(":")[0], "broken")
            hg38 = at
        elif st == "FAIL_hg38_to_hg19":
            cat, why = "not_liftable_back", bu.get(vid, (".", "unmapped"))[1]
            hg38 = f[5]
        elif st == "CHANGED_chrom":
            cat, why, hg38 = "chromosome_changed", "", f[5]
        else:
            cat, why, hg38 = "duplicated_region_position_changed", st, f[5]
        stats[o[0]][cat] += 1
        reasons[(cat, why.split(":")[0])] += 1
        u.write("\t".join(o + [cat, why, hg38 or ".", f[6] or "."]) + "\n")

ck = lambda c: (0, int(c[3:]), c) if c[3:].isdigit() else (1, {"X": 23, "Y": 24, "M": 25}.get(c[3:], 99), c)
matched.sort(key=lambda x: (ck(x[0]), x[1], x[2], x[3]))
with open(f"{prefix}_hg38.txt", "w") as m:
    for c, p, r, a, x in matched:
        m.write("\t".join([c, str(p), r, a] + x) + "\n")

with open(f"{prefix}_stats.tsv", "w") as s:
    s.write("chromosome\ttotal\tmatched\tmismatched\tpct_mismatched\t" + "\t".join(order) + "\n")
    tot = collections.Counter()
    for c in sorted(stats, key=ck):
        v = stats[c]; tot.update(v); mm = v["total"] - v["matched"]
        s.write(f"{c}\t{v['total']}\t{v['matched']}\t{mm}\t{100*mm/v['total']:.3f}\t" + "\t".join(str(v[k]) for k in order) + "\n")
    mm = tot["total"] - tot["matched"]
    s.write(f"ALL\t{tot['total']}\t{tot['matched']}\t{mm}\t{100*mm/tot['total']:.3f}\t" + "\t".join(str(tot[k]) for k in order) + "\n")
    s.write("\n#category\treason\tcount\n")
    for (cat, why), v in sorted(reasons.items(), key=lambda x: -x[1]):
        s.write(f"#{cat}\t{why or '-'}\t{v}\n")
print(open(f"{prefix}_stats.tsv").read())
