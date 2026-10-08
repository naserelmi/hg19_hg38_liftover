#!/usr/bin/env python3
"""Comparisons for the Iranome hg19 -> hg38 -> hg19 liftover check.

compare.py roundtrip <orig_norm.vcf> <back_norm.vcf> <fwd_unmap> <back_unmap> <out_prefix>
    Join by ID (original hg19 row) and classify every original variant.
compare.py provided <mine38_norm.vcf> <provided38_norm.vcf> <out_prefix>
    Compare our hg38 liftover against the provided hg38 table (key = chr,pos,ref,alt).
compare.py keys <a_norm.vcf> <b_norm.vcf> <out_prefix>
    Key-based (chr,pos,ref,alt) set comparison (provided hg38 lifted back vs original hg19).
"""
import sys, gzip, collections


def opener(p):
    return gzip.open(p, "rt") if p.endswith("gz") else open(p)


def read_vcf(p):
    """yield (chrom,pos,id,ref,alt,info-dict)"""
    with opener(p) as fh:
        for l in fh:
            if l[0] == "#":
                continue
            f = l.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "1") for kv in f[7].split(";"))
            yield f[0], int(f[1]), f[2], f[3], f[4], info


def vtype(ref, alt):
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(ref) == 1 and alt[0] == ref[0]:
        return "INS"
    if len(alt) == 1 and ref[0] == alt[0]:
        return "DEL"
    return "COMPLEX"


def read_unmap(p):
    """CrossMap .unmap: VCF lines with a reason in the last column."""
    d = {}
    try:
        with opener(p) as fh:
            for l in fh:
                if l[0] == "#":
                    continue
                f = l.rstrip("\n").split("\t")
                d[f[2]] = f[-1] if len(f) > 8 else "unmapped"
    except FileNotFoundError:
        pass
    return d


def roundtrip(orig, back, fwd_unmap, back_unmap, out, fwd=None):
    o = {r[2]: r for r in read_vcf(orig)}
    h38 = {r[2]: f"{r[0]}:{r[1]}:{r[3]}:{r[4]}" for r in read_vcf(fwd)} if fwd else {}
    b = collections.defaultdict(list)
    for r in read_vcf(back):
        b[r[2]].append(r)
    fu, bu = read_unmap(fwd_unmap), read_unmap(back_unmap)
    cnt = collections.Counter()
    bytype = collections.defaultdict(collections.Counter)
    with open(out + ".roundtrip_all.tsv", "w") as w, open(out + ".roundtrip_changed.tsv", "w") as wc:
        hdr = "ID\tstatus\ttype\thg19_table_row\thg19_normalized\thg38_lifted\thg19_roundtrip\tAF\tnote\n"
        w.write(hdr); wc.write(hdr)
        for vid, (c, p, _, ref, alt, info) in o.items():
            t = vtype(ref, alt)
            src = f"{c}:{p}:{ref}:{alt}"
            l38 = h38.get(vid, "")
            raw = info.get("ORIG", "")
            if vid in fu:
                st, dst, note = "FAIL_hg19_to_hg38", "", fu[vid]
            elif vid not in b:
                st, dst, note = "FAIL_hg38_to_hg19", "", bu.get(vid, "unmapped")
            else:
                rs = b[vid]
                dsts = [f"{x[0]}:{x[1]}:{x[3]}:{x[4]}" for x in rs]
                dst = ",".join(dsts)
                note = ""
                if src in dsts and len(dsts) == 1:
                    st = "IDENTICAL"
                elif src in dsts:
                    st = "IDENTICAL_plus_extra"
                else:
                    x = rs[0]
                    if x[0] != c:
                        st = "CHANGED_chrom"
                    elif x[1] != p and (x[3], x[4]) == (ref, alt):
                        st = "CHANGED_position"
                    elif x[1] == p and (x[3], x[4]) != (ref, alt):
                        st = "CHANGED_allele"
                    else:
                        st = "CHANGED_position_and_allele"
            cnt[st] += 1
            bytype[t][st] += 1
            line = f"{vid}\t{st}\t{t}\t{raw}\t{src}\t{l38}\t{dst}\t{info.get('AF','')}\t{note}\n"
            w.write(line)
            if st != "IDENTICAL":
                wc.write(line)
    print("== round trip (per original hg19 variant) ==")
    tot = sum(cnt.values())
    for k, v in cnt.most_common():
        print(f"{k:32s}{v:>10d}  {100*v/tot:6.2f}%")
    print(f"{'TOTAL':32s}{tot:>10d}")
    print("\nby variant type:")
    sts = [k for k, _ in cnt.most_common()]
    print("type\t" + "\t".join(sts))
    for t, c in sorted(bytype.items()):
        print(t + "\t" + "\t".join(str(c[s]) for s in sts))


def keyset(p):
    d = {}
    for c, pos, vid, ref, alt, info in read_vcf(p):
        d.setdefault((c, pos, ref, alt), (vid, info))
    return d


def provided(mine, prov, out):
    m, pv = keyset(mine), keyset(prov)
    # index provided by chrom,pos and by chrom,pos,AF to explain non-matches
    by_pos = collections.defaultdict(list)
    for k, (vid, info) in pv.items():
        by_pos[(k[0], k[1])].append((k, info.get("AF")))
    both = m.keys() & pv.keys()
    af_eq = sum(1 for k in both if m[k][1].get("AF") == pv[k][1].get("AF"))
    print("== our hg38 liftover vs provided hg38 table ==")
    print(f"our lifted variants          {len(m):>10d}")
    print(f"provided hg38 variants       {len(pv):>10d}")
    print(f"exact match chr/pos/ref/alt  {len(both):>10d}  (AF identical in {af_eq})")
    cat = collections.Counter()
    pair = collections.Counter()
    with open(out + ".ours_not_in_provided.tsv", "w") as w:
        w.write("ID\thg19_orig\tour_hg38\tour_type\tAF\tcategory\tprovided_same_pos\n")
        for k in m.keys() - pv.keys():
            vid, info = m[k]
            same = by_pos.get((k[0], k[1]), [])
            same_af = [s for s in same if s[1] == info.get("AF")]
            if same_af:
                cg = "provided_same_pos_same_AF_diff_allele"
                pair[(vtype(k[2], k[3]), vtype(same_af[0][0][2], same_af[0][0][3]))] += 1
            elif same:
                cg = "provided_same_pos_other_AF"
            else:
                cg = "absent_at_pos"
            cat[cg] += 1
            sp = ",".join(f"{s[0][2]}>{s[0][3]}|AF={s[1]}" for s in same)
            w.write(f"{vid}\t{info.get('ORIG','')}\t{k[0]}:{k[1]}:{k[2]}:{k[3]}\t{vtype(k[2],k[3])}\t{info.get('AF')}\t{cg}\t{sp}\n")
    print(f"ours not in provided         {len(m.keys()-pv.keys()):>10d}")
    for k, v in cat.most_common():
        print(f"   {k:42s}{v:>9d}")
    print("   type change ours -> provided (same pos, same AF):")
    for (a, b), v in pair.most_common(15):
        print(f"      {a:8s}-> {b:8s}{v:>9d}")
    with open(out + ".provided_not_in_ours.tsv", "w") as w:
        w.write("provided_row\thg38\ttype\tAF\tREFOK\n")
        n = 0
        for k in pv.keys() - m.keys():
            vid, info = pv[k]
            n += 1
            w.write(f"{info.get('ORIG')}\t{k[0]}:{k[1]}:{k[2]}:{k[3]}\t{vtype(k[2],k[3])}\t{info.get('AF')}\t{info.get('REFOK')}\n")
    print(f"provided not in ours         {n:>10d}")


def keys(a, b, out):
    A, B = keyset(a), keyset(b)
    both = A.keys() & B.keys()
    print(f"A={len(A)} B={len(B)} shared={len(both)} onlyA={len(A.keys()-B.keys())} onlyB={len(B.keys()-A.keys())}")
    with open(out + ".onlyA.tsv", "w") as w:
        for k in sorted(A.keys() - B.keys()):
            w.write("\t".join(map(str, k)) + f"\t{A[k][1].get('AF')}\n")
    with open(out + ".onlyB.tsv", "w") as w:
        for k in sorted(B.keys() - A.keys()):
            w.write("\t".join(map(str, k)) + f"\t{B[k][1].get('ORIG')}\t{B[k][1].get('AF')}\n")


if __name__ == "__main__":
    {"roundtrip": roundtrip, "provided": provided, "keys": keys}[sys.argv[1]](*sys.argv[2:])
