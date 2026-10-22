#!/usr/bin/env bash
# Round-trip liftover (hg19 -> hg38 -> hg19) of the normalized Iranome VCF with transanno.
# Usage: bash code/transanno_roundtrip.sh <label> <hg19ToHg38 chain> <hg38ToHg19 chain> [input normalized hg19 VCF] [out base dir]
#   e.g. ucsc_chain    ref/hg19ToHg38.over.chain.gz            ref/hg38ToHg19.over.chain.gz
#        minimap2_chain transanno/chains/hg19ToHg38.transanno.chain.gz transanno/chains/hg38ToHg19.transanno.chain.gz
# Input:  default work_data/orig19.norm.vcf (from run_pipeline.sh). Output: <out base, default transanno>/<label>/
set -euo pipefail
D=/root/naser/liftover; C=$D/code; R=$D/ref
LABEL=$1; FWD=$(readlink -f $2); BACK=$(readlink -f $3)
IN=$(readlink -f ${4:-$D/work_data/orig19.norm.vcf}); BASE=$(readlink -m ${5:-$D/transanno})
W=$BASE/$LABEL; mkdir -p $W/tmp; cd $W
source /data1/HUVYA/opt/miniforge3/etc/profile.d/conda.sh; conda activate transanno
PRIM='^chr([0-9]+|X|Y|M)$'

norm() { bcftools sort -m 8G -Ou -T $W/tmp/ 2>/dev/null | bcftools norm -f "$1" -c w -Ov 2> norm_$2.log \
    | awk 'BEGIN{OFS="\t"} /^#/{print;next} {$4=toupper($4); $5=toupper($5); print}'; }
fail2unmap() { # transanno fail VCF -> unmap table (VCF columns + REASON last)
  awk -F'\t' 'BEGIN{OFS="\t"} /^#/{next} {r="unmapped"; if (match($8,/FAILED_REASON=[^;]+/)) r=substr($8,RSTART+14,RLENGTH-14); print $0, r}' "$1"; }
split_primary() { # $1 in.vcf, $2 primary out, appends non-primary records to $3 (unmap)
  awk -F'\t' -v p="$PRIM" -v u="$3" 'BEGIN{OFS="\t"} /^#/{print; next} $1~p{print; next} {print $0, "non_primary_contig:" $1 >> u}' "$1" > "$2"; }
strip() { awk -F'\t' 'BEGIN{OFS="\t"} /^##INFO=<ID=(ORIGINAL_|REF_CHANGED|MULTIMAP)/{next} /^#/{print;next}
  {gsub(/;(ORIGINAL_[A-Z]+=[^;]*|REF_CHANGED|MULTIMAP=[^;]*)/,"",$8); print}'; }

echo "### [$LABEL] 1. hg19 -> hg38"
transanno liftvcf -r $R/hg19.fa -q $R/hg38.fa -c $FWD -v $IN -o fwd38.raw.vcf -f fwd38.fail.vcf 2> fwd.log
tail -3 fwd.log
fail2unmap fwd38.fail.vcf > fwd38.unmap
split_primary fwd38.raw.vcf fwd38.prim.vcf fwd38.unmap
grep -v '^#' fwd38.prim.vcf | awk -F'\t' '$8~/REF_CHANGED/{print $3}' > fwd_ref_changed.ids
norm $R/hg38.fa fwd38 < fwd38.prim.vcf | strip > fwd38.norm.vcf

echo "### [$LABEL] 2. hg38 -> hg19"
transanno liftvcf -r $R/hg38.fa -q $R/hg19.fa -c $BACK -v fwd38.norm.vcf -o back19.raw.vcf -f back19.fail.vcf 2> back.log
tail -3 back.log
fail2unmap back19.fail.vcf > back19.unmap
split_primary back19.raw.vcf back19.prim.vcf back19.unmap
norm $R/hg19.fa back19 < back19.prim.vcf | strip > back19.norm.vcf

echo "### [$LABEL] 3. round trip"
python3 $C/compare.py roundtrip $IN back19.norm.vcf fwd38.unmap back19.unmap rt fwd38.norm.vcf > roundtrip_summary.txt
cat roundtrip_summary.txt
echo "forward REF/ALT swaps (REF_CHANGED): $(wc -l < fwd_ref_changed.ids)"
