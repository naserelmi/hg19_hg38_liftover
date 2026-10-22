#!/usr/bin/env bash
# Build DataIranome-format files and ANNOVAR filter databases (with .idx) for hg19 and hg38.
# Usage: bash code/make_annovar_db.sh [matched hg38 table] [out dir]
#        default: results_strict/1_matched_variants_hg38_strict.tsv -> annovar/   (run from /root/naser/liftover)
#
# Sources
#   hg38: strict set from code/strict_filter.py (all checks agree), or another matched table given as $1
#   hg19: work_data/orig19.norm.vcf            (DataIranome.txt, left-aligned/normalized on hg19)
# Both are de-duplicated (exact copies dropped; variants with conflicting rows keep the highest-AF row,
# see collapse_dups.py) so each variant appears once.
#
# Output in annovar/
#   DataIranome_hg19.txt, DataIranome_hg38.txt   same layout as DataIranome.txt: chr pos ref alt het hom AF,
#                                                no header, '_' = empty allele
#   hg19_iranome.txt(.idx), hg38_iranome.txt(.idx)  ANNOVAR filter databases:
#       Chr Start End Ref Alt Iranome_AF Iranome_Het Iranome_Hom Iranome_AN ; '-' = empty allele,
#       chromosome without 'chr'; deletion Start..End = deleted bases, insertion Start = End = base before it
#   conflicting_rows_hg19.tsv, conflicting_rows_hg38.tsv   variants that had several rows with different counts
set -euo pipefail
D=/root/naser/liftover
C=$D/code
ANNOVAR=/data1/HUVYA/opt/anaconda/annovar
BIN=1000
MATCHED=$(readlink -f ${1:-$D/results_strict/1_matched_variants_hg38_strict.tsv})
O=$(readlink -m ${2:-$D/annovar}); mkdir -p $O; cd $O

echo "## DataIranome-format tables"
python3 $C/vcf2table.py $D/work_data/orig19.norm.vcf | cut -f1-7 > hg19.tmp
cut -f1-7 $MATCHED > hg38.tmp
for b in hg19 hg38; do
  python3 $C/collapse_dups.py $b.tmp DataIranome_$b.txt conflicting_rows_$b.tsv
  rm $b.tmp
done

echo "## ANNOVAR databases + index (bin $BIN)"
for b in hg19 hg38; do
  { printf '#Chr\tStart\tEnd\tRef\tAlt\tIranome_AF\tIranome_Het\tIranome_Hom\tIranome_AN\n'
    awk -F'\t' 'BEGIN{OFS="\t"} {
        c=$1; sub(/^chr/,"",c); r=($3=="_")?"-":$3; a=($4=="_")?"-":$4
        e=(r=="-")?$2:$2+length(r)-1; an=int(($5+2*$6)/$7+0.5)
        print c,$2,e,r,a,$7,$5,$6,an }' DataIranome_$b.txt; } > ${b}_iranome.txt
  perl $ANNOVAR/index_annovar.pl ${b}_iranome.txt $BIN > ${b}_iranome.txt.idx
done
wc -l DataIranome_hg19.txt DataIranome_hg38.txt
