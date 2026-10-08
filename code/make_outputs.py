#!/usr/bin/env python3
"""Split the DataIranome hg19 -> hg38 -> hg19 liftover into result tables.

Usage: make_outputs.py <work dir> <out dir>
Reads from <work dir>: orig19.norm.vcf, fwd38.norm.vcf, fwd38.unmap, back19.unmap, rt.roundtrip_all.tsv
Writes to <out dir>:
  1_matched_variants_hg38.tsv          round trip gives back the identical hg19 variant
  2_altswapped_variants_hg38.tsv       hg38 reference = hg19 ALT; REF/ALT swapped and AF/genotypes recalculated
  3_deleted_variants_hg38.tsv          region deleted / partially deleted in hg38 (no hg38 position)
  4_chrom_changed_variants_hg38.tsv    lifted to hg38, but lifting back lands on another chromosome
  5_broken_variants_hg38.tsv           alt/unplaced contig only, split, span length changed, or not liftable back
  6_duplicate_regions_mismatch_hg38.tsv  lifted to hg38, but lifting back lands elsewhere on the same chromosome
  unmatched_variants_hg19.tsv          all of 2-6 in one list (input for realign_unmatched.sh)
  summary.txt
Coordinates/alleles use the DataIranome convention: '_' = empty allele, deletion Pos = first deleted base,
insertion Pos = base before the insertion. hg38 variants are left-aligned on hg38.
"""
import sys, os, collections

work, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]


def vcf(p):
    with open(p) as fh:
        for l in fh:
            if l[0] != "#":
                f = l.rstrip("\n").split("\t")
                yield f, dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)


def annovar(c, p, ref, alt):
    """VCF allele -> DataIranome style"""
    p = int(p)
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt:
        return c, p + 1, ref[1:], "_"
    if len(alt) > 1 and len(ref) == 1 and alt[0] == ref:
        return c, p, "_", alt[1:]
    return c, p, ref, alt


def hg19_row(info):
    c, s, e, r, a = info["ORIG"].split(":")
    return [f"chr{c}", s, r.replace("-", "_"), a.replace("-", "_")]


def unmap(p):
    d = {}
    with open(p) as fh:
        for l in fh:
            if l[0] == "#":
                continue
            f = l.rstrip("\n").split("\t")
            d[f[2]] = (f, f[-2], f[-1])           # record, LIFTED_AT, REASON
    return d


orig = {f[2]: (f, i) for f, i in vcf(f"{work}/orig19.norm.vcf")}
fwd = {f[2]: (f, i) for f, i in vcf(f"{work}/fwd38.norm.vcf")}
fu, bu = unmap(f"{work}/fwd38.unmap"), unmap(f"{work}/back19.unmap")
rt = {}
with open(f"{work}/rt.roundtrip_all.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        rt[f[0]] = f

# other ALT alleles at the same original hg19 site (for multi-allelic swap recalculation)
site = collections.defaultdict(list)
for vid, (f, i) in orig.items():
    c, s, e, r, a = i["ORIG"].split(":")
    site[(c, s, r)].append(vid)

G = ["Het", "Hom", "AF", "AC", "AN"]
H38 = ["Chr", "Pos", "Ref", "Alt"]
H19 = ["hg19_Chr", "hg19_Pos", "hg19_Ref", "hg19_Alt"]
files = {
    1: ("1_matched_variants_hg38.tsv", H38 + G + H19 + ["Strand", "hg38_ref_differs_from_hg19"]),
    2: ("2_altswapped_variants_hg38.tsv", H38 + G + H19 + ["Strand", "hg19_Het", "hg19_Hom", "hg19_AF", "hg19_AC", "multiallelic_site", "note"]),
    3: ("3_deleted_variants_hg38.tsv", H19 + G + ["reason"]),
    4: ("4_chrom_changed_variants_hg38.tsv", H38 + G + H19 + ["hg19_after_roundtrip"]),
    5: ("5_broken_variants_hg38.tsv", H38 + G + H19 + ["reason", "lifted_to"]),
    6: ("6_duplicate_regions_mismatch_hg38.tsv", H38 + G + H19 + ["hg19_after_roundtrip", "shift_bp"]),
}
fh = {k: open(f"{out}/{n}", "w") for k, (n, h) in files.items()}
for k, (n, h) in files.items():
    fh[k].write("#" + "\t".join(h) + "\n")
cnt = collections.Counter()
sub = collections.defaultdict(collections.Counter)
w = lambda k, row: fh[k].write("\t".join(map(str, row)) + "\n")

for vid, (f, i) in orig.items():
    g = [i["HET"], i["HOM"], i["AF"], i["AC"], i["AN"]]
    h19 = hg19_row(i)
    st = rt[vid][1]
    if st == "IDENTICAL":
        x, xi = fwd[vid]
        w(1, list(annovar(x[0], x[1], x[3], x[4])) + g + h19 + [xi.get("STRAND"), "yes" if "REFDIFF" in xi else "no"])
        cnt[1] += 1
    elif st == "FAIL_hg19_to_hg38":
        rec, at, why = fu[vid]
        if why == "ref_alt_swap":
            c, s, e, strand, newref = at.split(":")
            oldref = f[3] if strand == "+" else rc(f[3])
            het, hom, ac, an = int(i["HET"]), int(i["HOM"]), int(i["AC"]), int(i["AN"])
            c19, s19, _, r19, _ = i["ORIG"].split(":")
            others = [o for o in site[(c19, s19, r19)] if o != vid]
            if not others:
                nhet, nhom = het, an // 2 - het - hom      # hom-ref people become hom-alt
                nac = an - ac
                multi, note = "no", "Het unchanged; Hom = N_called - Het - Hom(hg19); AC = AN - AC; AF = 1 - AF"
            else:
                oac = sum(int(orig[o][1]["AC"]) for o in others)
                nac = max(an - ac - oac, 0)
                nhet, nhom = "NA", "NA"
                multi = "yes"
                note = f"AC = AN - AC - AC(other ALT alleles at site: {oac}); genotype counts not derivable"
            naf = f"{nac/an:.6g}" if an else "NA"
            w(2, list(annovar(c, s, newref, oldref)) + [nhet, nhom, naf, nac, an] + h19 +
              [strand, het, hom, i["AF"], ac, multi, note])
            cnt[2] += 1; sub[2][multi] += 1
        elif why in ("Deleted_in_new", "Partially_deleted_in_new"):
            w(3, h19 + g + [why]); cnt[3] += 1; sub[3][why] += 1
        else:
            w(5, [".", ".", ".", "."] + g + h19 + [why, at]); cnt[5] += 1; sub[5][why.split(":")[0]] += 1
    elif st == "FAIL_hg38_to_hg19":
        x, xi = fwd[vid]
        w(5, list(annovar(x[0], x[1], x[3], x[4])) + g + h19 + ["hg38_not_liftable_back:" + bu.get(vid, (0, 0, "unmapped"))[2], "."])
        cnt[5] += 1; sub[5]["hg38_not_liftable_back"] += 1
    elif st == "CHANGED_chrom":
        x, xi = fwd[vid]
        w(4, list(annovar(x[0], x[1], x[3], x[4])) + g + h19 + [rt[vid][6]]); cnt[4] += 1
    else:   # CHANGED_position / CHANGED_position_and_allele
        x, xi = fwd[vid]
        back = rt[vid][6].split(",")[0].split(":")
        shift = int(back[1]) - int(f[1])
        w(6, list(annovar(x[0], x[1], x[3], x[4])) + g + h19 + [rt[vid][6], shift]); cnt[6] += 1; sub[6][st] += 1

for f_ in fh.values():
    f_.close()
tot = sum(cnt.values())
with open(f"{out}/summary.txt", "w") as s:
    s.write(f"DataIranome hg19 variants: {tot}\n")
    for k, (n, h) in files.items():
        s.write(f"{n:42s}{cnt[k]:>10d}  {100*cnt[k]/tot:6.3f}%   {dict(sub[k]) if sub[k] else ''}\n")
print(open(f"{out}/summary.txt").read())

# combined list of every variant that did not round-trip identically (input for realign_unmatched.sh)
with open(f"{out}/unmatched_variants_hg19.tsv", "w") as u:
    u.write("#hg19_Chr\thg19_Pos\thg19_Ref\thg19_Alt\tcategory\treason\thg38_Chr\thg38_Pos\n")
    for k in (2, 3, 4, 5, 6):
        name = files[k][0]
        cat_ = name.split("_", 1)[1].replace("_hg38.tsv", "")
        with open(f"{out}/{name}") as fh2:
            hdr = fh2.readline().lstrip("#").rstrip("\n").split("\t")
            ix = {h: n for n, h in enumerate(hdr)}
            for l in fh2:
                f = l.rstrip("\n").split("\t")
                why = f[ix["reason"]] if "reason" in ix else cat_
                c38, p38 = (f[ix["Chr"]], f[ix["Pos"]]) if "Chr" in ix else (".", ".")
                u.write("\t".join([f[ix["hg19_Chr"]], f[ix["hg19_Pos"]], f[ix["hg19_Ref"]], f[ix["hg19_Alt"]], cat_, why, c38, p38]) + "\n")
