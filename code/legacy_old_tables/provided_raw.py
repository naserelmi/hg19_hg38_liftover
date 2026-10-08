#!/usr/bin/env python3
"""Row-level audit of the provided hg38 table against the hg19 table.

Each hg19 row's Start..End is lifted with UCSC liftOver; the provided hg38 row at the lifted
Start (same chr) with the same AF is taken as its counterpart, and the Ref/Alt pair is compared.
Usage: provided_raw.py <hg19 table> <hg38 table> <hg19ToHg38 chain> <out prefix>
"""
import sys, subprocess, collections, os

t19, t38, chain, out = sys.argv[1:5]


def rows(p):
    with open(p) as fh:
        for n, l in enumerate(fh, 1):
            if l[0] != "#":
                c, s, e, r, a, af = l.rstrip("\n").split("\t")[:6]
                yield n, c, int(s), int(e), r, a, af


def kind(r, a):
    if a == "-": return "DEL"
    if r == "-": return "INS"
    if len(r) == len(a) == 1: return "SNV"
    if len(r) == 1 and a.startswith(r): return "INS_anchored"
    return "OTHER"


p38 = collections.defaultdict(list)          # (chr,start) -> rows
p38_all = {}
for n, c, s, e, r, a, af in rows(t38):
    p38[(c, s)].append((n, r, a, af))

bed, lo, un = out + ".in.bed", out + ".out.bed", out + ".un.bed"
h19 = {}
with open(bed, "w") as b:
    for n, c, s, e, r, a, af in rows(t19):
        h19[n] = (c, s, e, r, a, af)
        b.write(f"chr{c}\t{s-1}\t{max(e, s)}\t{n}\n")
subprocess.run(["liftOver", "-minMatch=1", bed, chain, lo, un], check=True, stderr=subprocess.DEVNULL)
lifted = {}
with open(lo) as fh:
    for l in fh:
        c, s, e, n = l.split("\t")[:4]
        lifted[int(n)] = (c[3:], int(s) + 1)

cat = collections.Counter()
trans = collections.Counter()
used = set()
ex = collections.defaultdict(list)
with open(out + ".rows.tsv", "w") as w:
    w.write("hg19_row\thg19\thg38_lifted_start\tprovided_hg38\tcategory\n")
    for n, (c, s, e, r, a, af) in h19.items():
        k19 = kind(r, a)
        if n not in lifted:
            cg, prov = "hg19_row_not_liftable", ""
        else:
            lc, ls = lifted[n]
            cand = [x for x in p38.get((lc, ls), []) if x[3] == af]
            if not cand:
                # was it copied with hg19 coordinates instead of lifted ones?
                same19 = [x for x in p38.get((c, s), []) if x[3] == af] if (c, s) != (lc, ls) else []
                cg = "missing_in_provided" if not same19 else "provided_at_hg19_coordinate"
                prov = ";".join(f"{x[1]}>{x[2]}" for x in same19)
            else:
                exact = [x for x in cand if (x[1], x[2]) == (r, a)]
                x = exact[0] if exact else cand[0]
                used.add(x[0])
                prov = f"{lc}:{ls}:{x[1]}>{x[2]}"
                if exact:
                    cg = "identical"
                else:
                    cg = f"changed:{k19}->{kind(x[1], x[2])}"
                    # sub-classify
                    if k19 == "SNV" and x[2] == "-" and x[1] == r:
                        cg += " (ALT allele lost, SNV became 1-bp deletion)"
                    elif k19 == "DEL" and x[2] == "-" and r.startswith(x[1]):
                        cg += " (deleted sequence truncated to first base)"
                    elif k19 == "INS" and x[1] != "-" and x[2] == x[1] + a[1:] and a[0] == x[1]:
                        cg += " (anchor base merged into inserted sequence: 1 base lost)"
                    elif k19 == "INS" and x[1] != "-" and x[2] == x[1] + a:
                        cg += " (anchored correctly)"
            ls_ = lifted[n][1] if n in lifted else ""
        cat[cg] += 1
        if len(ex[cg]) < 5:
            ex[cg].append(f"{c}:{s}-{e} {r}>{a} AF={af}  ->  {prov}")
        w.write(f"{n}\t{c}:{s}:{e}:{r}:{a}:{af}\t{ls_}\t{prov}\t{cg}\n")

orphan = [n for k, v in p38.items() for (n, r, a, af) in v if n not in used]
print(f"hg19 rows {len(h19)}, provided hg38 rows {sum(len(v) for v in p38.values())}")
for k, v in sorted(cat.items(), key=lambda x: -x[1]):
    print(f"{v:>9d}  {k}")
    for e in ex[k][:3]:
        print("              e.g.", e)
print(f"provided hg38 rows with no hg19 counterpart: {len(orphan)}")
for f in (bed, lo, un):
    os.remove(f)
