#!/usr/bin/env bash
# DataIranome (hg19) -> hg38 -> hg19 liftover with round-trip validation.
# Usage: bash code/run_pipeline.sh            (run from /root/naser/liftover)
# Needs: ref/{hg19,hg38}.fa.gz (bgzip + faidx), ref/hg19ToHg38.over.chain.gz, ref/hg38ToHg19.over.chain.gz,
#        conda env `liftover` (python3 + pysam + ucsc-liftover), bcftools >= 1.13.
set -euo pipefail
D=/root/naser/liftover
C=$D/code; R=$D/ref; W=$D/work_data; O=$D/results
export PATH=/data1/HUVYA/opt/miniforge3/envs/liftover/bin:$PATH
mkdir -p $W/tmp $O; cd $W

addtag() { # copy CHROM:POS:REF:ALT into INFO/<tag>
  awk -v t="$1" 'BEGIN{OFS="\t"} /^##INFO/&&!d{print "##INFO=<ID=" t ",Number=1,Type=String,Description=\"position before this step\">";d=1} /^#/{print;next} {$8=$8";"t"="$1":"$2":"$4":"$5; print}'
}
norm() { # $1 genome, $2 log name: sort, left-align, normalize, upper-case
  bcftools sort -Ou -T $W/tmp/ 2>/dev/null | bcftools norm -f "$1" -c w -Ov 2> norm_$2.log \
    | awk 'BEGIN{OFS="\t"} /^#/{print;next} {$4=toupper($4); $5=toupper($5); print}'
  grep -E '^Lines' norm_$2.log >&2 || true
}

echo "### 0. DataIranome -> 6-column table (+HET,HOM,AC,AN); junk rows (non-numeric counts) dropped"
# DataIranome columns: chr pos ref alt het hom AF ; '_' = empty allele ; AN recovered as (het+2hom)/AF
awk -F'\t' 'BEGIN{OFS="\t"; print "#Chr","Start","End","Ref","Alt","AF","HET","HOM","AC","AN"}
  $5~/^[0-9]+$/ { c=$1; sub(/^chr/,"",c); r=($3=="_")?"-":$3; a=($4=="_")?"-":$4;
     e=(r=="-")?$2:$2+length(r)-1; ac=$5+2*$6; an=int(ac/$7+0.5); print c,$2,e,r,a,$7,$5,$6,ac,an }' \
  $D/DataIranome.txt > data19.tsv
echo "rows: $(($(wc -l < data19.tsv)-1))"

echo "### 1. table -> VCF (anchor bases from hg19, REF checked), normalize on hg19"
python3 $C/annovar2vcf.py data19.tsv $R/hg19.fa.gz orig19.raw.vcf d19
addtag RAW19 < orig19.raw.vcf | norm $R/hg19.fa.gz orig19 > orig19.norm.vcf

echo "### 2. hg19 -> hg38 (UCSC liftOver + chain, REF re-read from hg38), normalize on hg38"
python3 $C/liftvcf.py orig19.norm.vcf $R/hg19ToHg38.over.chain.gz $R/hg38.fa.gz fwd38.vcf fwd38.unmap
norm $R/hg38.fa.gz fwd38 < fwd38.vcf | addtag HG38 > fwd38.norm.vcf

echo "### 3. hg38 -> hg19, normalize on hg19"
python3 $C/liftvcf.py fwd38.norm.vcf $R/hg38ToHg19.over.chain.gz $R/hg19.fa.gz back19.vcf back19.unmap
norm $R/hg19.fa.gz back19 < back19.vcf > back19.norm.vcf

echo "### 4. round-trip classification (per original variant)"
python3 $C/compare.py roundtrip orig19.norm.vcf back19.norm.vcf fwd38.unmap back19.unmap rt fwd38.norm.vcf | tee roundtrip_summary.txt

echo "### 5. result tables"
python3 $C/make_outputs.py $W $O
python3 $C/vcf2table.py fwd38.norm.vcf > $O/all_lifted_variants_hg38.tsv
bgzip -c fwd38.norm.vcf > $O/all_lifted_variants_hg38.vcf.gz && tabix -f -p vcf $O/all_lifted_variants_hg38.vcf.gz
gzip -c rt.roundtrip_all.tsv > $O/roundtrip_per_variant.tsv.gz
cp roundtrip_summary.txt $O/
