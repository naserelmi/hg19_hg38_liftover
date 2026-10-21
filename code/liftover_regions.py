#!/usr/bin/env python3
"""Round-trip liftover (hg19 -> hg38 -> hg19) of a region table (BED-like, with a header).

Usage:
  liftover_regions.py IN_TABLE OUT_PREFIX [--chr chr_num --start start --end end]
                      [--one-based] [--chain-fwd F] [--chain-back B] [--min-match 1]

A region is MATCHED when it lifts to a primary hg38 chromosome and lifts back to exactly the same
chr/start/end. Coordinates are 0-based half-open (BED/UCSC) unless --one-based.

Output
  OUT_PREFIX_hg38.txt       matched regions only, same header/columns/order as the input; the three
                            coordinate columns hold the hg38 position, every other column is copied as is
  OUT_PREFIX_unmatched.tsv  every other region: original columns + status, reason, hg38 region,
                            hg19 region after the round trip
  OUT_PREFIX_stats.tsv      per hg19 chromosome: total, matched, mismatched, % mismatched, and why
"""
import argparse, collections, os, re, subprocess, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("table"); ap.add_argument("prefix")
ap.add_argument("--chr", default="chr_num"); ap.add_argument("--start", default="start"); ap.add_argument("--end", default="end")
ap.add_argument("--one-based", action="store_true")
ap.add_argument("--chain-fwd", default="/root/naser/liftover/ref/hg19ToHg38.over.chain.gz")
ap.add_argument("--chain-back", default="/root/naser/liftover/ref/hg38ToHg19.over.chain.gz")
ap.add_argument("--min-match", default="1")
a = ap.parse_args()
PRIMARY = re.compile(r"^chr([0-9]+|X|Y|M)$")

with open(a.table) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
ci, si, ei = header.index(a.chr), header.index(a.start), header.index(a.end)
off = 1 if a.one_based else 0


def lift(bed_rows, chain, tmp, tag):
    """bed_rows: (id, chr, start0, end) -> ({id: (chr, start0, end)}, {id: reason})"""
    i, o, u = (os.path.join(tmp, f"{tag}.{x}") for x in ("in", "out", "un"))
    with open(i, "w") as f:
        for rid, c, s, e in bed_rows:
            f.write(f"{c}\t{s}\t{e}\t{rid}\n")
    subprocess.run(["liftOver", f"-minMatch={a.min_match}", i, chain, o, u], check=True, stderr=subprocess.DEVNULL)
    ok = {}
    for l in open(o):
        c, s, e, rid = l.rstrip("\n").split("\t")[:4]
        ok[int(rid)] = (c, int(s), int(e))
    why, reason = "unmapped", {}
    for l in open(u):
        if l.startswith("#"):
            why = l[1:].strip().replace(" ", "_")
        else:
            reason[int(l.split("\t")[3])] = why
    return ok, reason


src = [(n, r[ci], int(r[si]) - off, int(r[ei])) for n, r in enumerate(rows)]
with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.prefix)) or ".") as tmp:
    fwd, fwd_un = lift(src, a.chain_fwd, tmp, "fwd")
    prim = {n: v for n, v in fwd.items() if PRIMARY.match(v[0])}
    back, back_un = lift([(n, *v) for n, v in prim.items()], a.chain_back, tmp, "back")

fmt = lambda v: f"{v[0]}:{v[1] + off}-{v[2]}" if v else "."
stats = collections.defaultdict(collections.Counter)
reasons = []
lendiff = 0
with open(f"{a.prefix}_hg38.txt", "w") as m, open(f"{a.prefix}_unmatched.tsv", "w") as u:
    m.write("\t".join(header) + "\n")
    u.write("\t".join(header + ["liftover_status", "reason", "hg38_region", "hg19_after_roundtrip"]) + "\n")
    for n, c, s, e in src:
        r = rows[n]
        if n in fwd_un:
            st, why = "not_liftable_to_hg38", fwd_un[n]
        elif n not in prim:
            st, why = "hg38_non_primary_contig", fwd[n][0]
        elif n in back_un:
            st, why = "not_liftable_back_to_hg19", back_un[n]
        elif back[n] != (c, s, e):
            b = back[n]
            st = "roundtrip_changed"
            why = "other_chromosome" if b[0] != c else ("same_start_other_end" if b[1] == s else
                  "same_end_other_start" if b[2] == e else "shifted")
        else:
            st, why = "matched", ""
        stats[c]["total"] += 1
        stats[c][st] += 1
        if st == "matched":
            h = fwd[n]
            lendiff += (h[2] - h[1]) != (e - s)
            out = list(r)
            out[ci], out[si], out[ei] = h[0], str(h[1] + off), str(h[2])
            m.write("\t".join(out) + "\n")
        else:
            u.write("\t".join(r + [st, why, fmt(fwd.get(n)), fmt(back.get(n))]) + "\n")
            reasons.append((st, why))

cats = ["not_liftable_to_hg38", "hg38_non_primary_contig", "not_liftable_back_to_hg19", "roundtrip_changed"]
ck = lambda c: (0, int(c[3:]), c) if c[3:].isdigit() else (1, {"X": 23, "Y": 24, "M": 25}.get(c[3:], 99), c)
with open(f"{a.prefix}_stats.tsv", "w") as f:
    f.write("hg19_chromosome\ttotal\tmatched\tmismatched\tpct_mismatched\t" + "\t".join(cats) + "\n")
    tot, ptot = collections.Counter(), collections.Counter()
    for c in sorted(stats, key=ck):
        s = stats[c]
        mm = s["total"] - s["matched"]
        f.write(f"{c}\t{s['total']}\t{s['matched']}\t{mm}\t{100*mm/s['total']:.2f}\t" + "\t".join(str(s[k]) for k in cats) + "\n")
        tot.update(s)
        if PRIMARY.match(c):
            ptot.update(s)
    for name, s in (("ALL_primary_chromosomes", ptot), ("ALL", tot)):
        mm = s["total"] - s["matched"]
        f.write(f"{name}\t{s['total']}\t{s['matched']}\t{mm}\t{100*mm/s['total']:.2f}\t" + "\t".join(str(s[k]) for k in cats) + "\n")

print(open(f"{a.prefix}_stats.tsv").read())
print("mismatch reasons:")
for (st, why), v in collections.Counter((s, w.split(":")[0] if s != "hg38_non_primary_contig" else "non_primary")
                                         for s, w in reasons).most_common():
    print(f"  {v:>7d}  {st}  {why}")
print(f"matched regions whose hg38 length differs from hg19: {lendiff}")
