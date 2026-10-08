#!/usr/bin/env python3
"""Normalized VCF -> DataIranome-style table (Chr Pos Ref Alt Het Hom AF AC AN), '_' for empty allele.
Indels are written ANNOVAR style: deletion Pos = first deleted base, insertion Pos = base before the insertion."""
import sys
print("#Chr\tPos\tRef\tAlt\tHet\tHom\tAF\tAC\tAN\thg19_source")
for l in open(sys.argv[1]):
    if l[0] == "#":
        continue
    c, p, vid, ref, alt, _, _, info = l.rstrip("\n").split("\t")[:8]
    d = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
    p = int(p)
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt:
        p, ref, alt = p + 1, ref[1:], "_"
    elif len(alt) > 1 and len(ref) == 1 and alt[0] == ref:
        ref, alt = "_", alt[1:]
    print(f"{c}\t{p}\t{ref}\t{alt}\t{d['HET']}\t{d['HOM']}\t{d['AF']}\t{d['AC']}\t{d['AN']}\t{d['ORIG']}")
