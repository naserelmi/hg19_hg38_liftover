#!/usr/bin/env python3
"""Collapse duplicate variants in a DataIranome-style table (chr pos ref alt het hom AF ...).

Usage: collapse_dups.py <in.tsv> <out.txt> <conflicts.tsv>
 - exact duplicate rows are dropped;
 - when one variant (chr,pos,ref,alt) still has several rows with different counts (the same indel
   written at different positions in a repeat, merged by left-alignment), the row with the highest AF
   is kept and all of its rows are written to <conflicts.tsv>.
Writes 7 columns (chr pos ref alt het hom AF), no header, sorted by chromosome and position.
"""
import sys, re, collections

src, out, conf = sys.argv[1:4]
rows = collections.defaultdict(set)
for l in open(src):
    if l[0] == "#":
        continue
    f = l.rstrip("\n").split("\t")
    rows[tuple(f[:4])].add(tuple(f[4:7]))

ck = lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, {"X": 23, "Y": 24, "M": 25}.get(c[3:], 99), c)
keys = sorted(rows, key=lambda k: (ck(k[0]), int(k[1]), k[2], k[3]))
n_conf = 0
with open(out, "w") as o, open(conf, "w") as c:
    c.write("#Chr\tPos\tRef\tAlt\tHet\tHom\tAF\tkept\n")
    for k in keys:
        v = rows[k]
        best = max(v, key=lambda x: float(x[2]))
        if len(v) > 1:
            n_conf += 1
            for x in sorted(v, key=lambda x: -float(x[2])):
                c.write("\t".join(k + x) + ("\tyes\n" if x == best else "\tno\n"))
        o.write("\t".join(k + best) + "\n")
print(f"{src}: distinct variants {len(keys)}, collapsed conflicting variants {n_conf}", file=sys.stderr)
