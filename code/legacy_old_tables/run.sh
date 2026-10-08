#!/usr/bin/env bash
# Iranome liftover check: hg19 -> hg38 -> hg19 round trip, and comparison with the provided hg38 table.
set -euo pipefail
D=/root/naser/liftover
R=$D/ref
W=$D/work
export PATH=/data1/HUVYA/opt/miniforge3/envs/liftover/bin:$PATH
mkdir -p $W/tmp; cd $W; cp $D/*.py .

addtag() { # copy CHROM:POS:REF:ALT into INFO/<tag>
  awk -v t="$1" 'BEGIN{OFS="\t"} /^##INFO/&&!d{print "##INFO=<ID=" t ",Number=1,Type=String,Description=\"position before this step\">";d=1} /^#/{print;next} {$8=$8";"t"="$1":"$2":"$4":"$5; print}'
}
norm() { # $1 genome, $2 log name; stdin vcf -> stdout sorted, left-aligned, normalized, upper-case vcf
  bcftools sort -Ou -T $W/tmp/ 2>/dev/null | bcftools norm -f "$1" -c w -Ov 2> norm_$2.log \
    | awk 'BEGIN{OFS="\t"} /^#/{print;next} {$4=toupper($4); $5=toupper($5); print}'
  grep -E '^Lines' norm_$2.log >&2 || true
}

echo "### 1. hg19 table -> VCF, left-align/normalize on hg19"
python3 annovar2vcf.py $D/hg19_iranome.txt $R/hg19.fa.gz orig19.raw.vcf h19
addtag RAW19 < orig19.raw.vcf | norm $R/hg19.fa.gz orig19 > orig19.norm.vcf

echo "### 2. hg19 -> hg38 (UCSC liftOver + chain), re-normalize on hg38"
python3 liftvcf.py orig19.norm.vcf $R/hg19ToHg38.over.chain.gz $R/hg38.fa.gz fwd38.vcf fwd38.unmap
norm $R/hg38.fa.gz fwd38 < fwd38.vcf | addtag HG38 > fwd38.norm.vcf

echo "### 3. hg38 -> hg19 (UCSC liftOver + chain), re-normalize on hg19"
python3 liftvcf.py fwd38.norm.vcf $R/hg38ToHg19.over.chain.gz $R/hg19.fa.gz back19.vcf back19.unmap
norm $R/hg19.fa.gz back19 < back19.vcf > back19.norm.vcf

echo "### 4. round-trip classification"
python3 compare.py roundtrip orig19.norm.vcf back19.norm.vcf fwd38.unmap back19.unmap rt fwd38.norm.vcf | tee roundtrip_summary.txt

echo "### 5. provided hg38 table -> VCF, normalize, compare with our liftover"
python3 annovar2vcf.py $D/hg38_iranome.txt $R/hg38.fa.gz prov38.raw.vcf p38
norm $R/hg38.fa.gz prov38 < prov38.raw.vcf > prov38.norm.vcf
python3 compare.py provided fwd38.norm.vcf prov38.norm.vcf prov | tee provided_summary.txt

echo "### 6. provided hg38 -> hg19, compare with original hg19 (key based)"
python3 liftvcf.py prov38.norm.vcf $R/hg38ToHg19.over.chain.gz $R/hg19.fa.gz prov_back19.vcf prov_back19.unmap
norm $R/hg19.fa.gz provback < prov_back19.vcf > prov_back19.norm.vcf
python3 compare.py keys orig19.norm.vcf prov_back19.norm.vcf provback | tee provback_summary.txt
