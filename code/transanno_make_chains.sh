#!/usr/bin/env bash
# transanno chains from minimap2 alignments (README: minimap2 -cx asm5 --cs QUERY REFERENCE; transanno minimap2chain)
set -euo pipefail
source /data1/HUVYA/opt/miniforge3/etc/profile.d/conda.sh; conda activate transanno
cd /root/naser/liftover/ref
for g in hg19 hg38; do
  [ -s ${g}.fa ] || { zcat ${g}.fa.gz > ${g}.fa; samtools faidx ${g}.fa; }
  [ -s ${g}_primary.fa ] || { samtools faidx ${g}.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > ${g}_primary.fa; samtools faidx ${g}_primary.fa; }
done
cd ../transanno/chains
( minimap2 -t 90 -cx asm5 --cs ../../ref/hg38_primary.fa ../../ref/hg19_primary.fa > hg19ToHg38.paf 2> mm_19to38.log && transanno minimap2chain hg19ToHg38.paf --output hg19ToHg38.transanno.chain && gzip -f hg19ToHg38.transanno.chain && echo FWD_DONE ) &
( minimap2 -t 90 -cx asm5 --cs ../../ref/hg19_primary.fa ../../ref/hg38_primary.fa > hg38ToHg19.paf 2> mm_38to19.log && transanno minimap2chain hg38ToHg19.paf --output hg38ToHg19.transanno.chain && gzip -f hg38ToHg19.transanno.chain && echo BACK_DONE ) &
wait; echo ALL_DONE
