#!/usr/bin/env bash
# Round-trip liftover (hg19 -> hg38 -> hg19) of a VCF-style variant table without header:
#   chr pos ref alt [extra columns...]        (e.g. AllVcf_V2.txt)
# Usage: bash code/vcftab_roundtrip.sh <table> <out dir>
# Output (<out dir>/<name>_*): _hg38.txt (exact matches, same structure, hg38), _unmatched.tsv, _stats.tsv
set -euo pipefail
D=/root/naser/liftover; C=$D/code; R=$D/ref
IN=$(readlink -f "$1"); OUT=$(readlink -m "$2"); NAME=$(basename "${IN%.*}")
W=$OUT/work; mkdir -p $W/tmp
export PATH=/data1/HUVYA/opt/miniforge3/envs/liftover/bin:$PATH
cd $W

addtag() { awk -v t="$1" 'BEGIN{OFS="\t"} /^##INFO/&&!d{print "##INFO=<ID=" t ",Number=1,Type=String,Description=\"position before this step\">";d=1} /^#/{print;next} {$8=$8";"t"="$1":"$2":"$4":"$5; print}'; }
norm() { bcftools sort -m 8G -Ou -T $W/tmp/ 2>/dev/null | bcftools norm -f "$1" -c w -Ov 2> norm_$2.log \
    | awk 'BEGIN{OFS="\t"} /^#/{print;next} {$4=toupper($4); $5=toupper($5); print}'; grep -E '^Lines' norm_$2.log >&2 || true; }

echo "### 1. table -> VCF (REF checked against hg19), normalize on hg19"
python3 $C/vcftab2vcf.py $IN $R/hg19.fa.gz orig19.raw.vcf
addtag RAW19 < orig19.raw.vcf | norm $R/hg19.fa.gz orig19 > orig19.norm.vcf
echo "### 2. hg19 -> hg38"
python3 $C/liftvcf.py orig19.norm.vcf $R/hg19ToHg38.over.chain.gz $R/hg38.fa.gz fwd38.vcf fwd38.unmap
norm $R/hg38.fa.gz fwd38 < fwd38.vcf | addtag HG38 > fwd38.norm.vcf
echo "### 3. hg38 -> hg19"
python3 $C/liftvcf.py fwd38.norm.vcf $R/hg38ToHg19.over.chain.gz $R/hg19.fa.gz back19.vcf back19.unmap
norm $R/hg19.fa.gz back19 < back19.vcf > back19.norm.vcf
echo "### 4. round trip"
python3 $C/compare.py roundtrip orig19.norm.vcf back19.norm.vcf fwd38.unmap back19.unmap rt fwd38.norm.vcf > roundtrip_summary.txt
cat roundtrip_summary.txt
echo "### 5. outputs"
python3 $C/vcftab_outputs.py $W $IN $OUT/$NAME
