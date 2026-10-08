#!/usr/bin/env bash
# realign_unmatched.sh - pull the reads behind hg19 variants that did not lift over cleanly
# out of an hg19 BAM/CRAM, write them as paired FASTQ, and realign them to hg38.
#
# Usage:
#   realign_unmatched.sh -v VARIANTS [-v VARIANTS2 ...] -i SAMPLE.cram|bam [-i ...] [options]
#
#   -v FILE   hg19 variant list. Columns are found by header name: hg19_Chr/hg19_Pos/hg19_Ref
#             (the results/*.tsv tables and results/unmatched_variants_hg19.tsv), or Chr/Pos/Ref.
#             A file without a header is read as chr, pos, ref, alt. '_' = empty allele.
#   -i FILE   hg19-aligned BAM or CRAM (indexed). Repeat -i or pass a text file with one path per line.
#   -o DIR    output directory                                  [./realigned_hg38]
#   -p INT    padding around each variant, bp                   [500]
#   -t INT    threads                                           [32]
#   -r FASTA  hg19 reference (CRAM decoding)                    [/root/naser/liftover/ref/hg19.fa.gz]
#   -x PREFIX hg38 bwa index prefix                             [/root/naser/liftover/ref/hg38_noalt_bwa/hg38_noalt.fa]
#   -f        fast mode: take only reads inside the regions (no whole-file scan for their mates;
#             a pair whose mate lies outside becomes a singleton)
#   -k        keep intermediate files
#
# Output per sample (<name> = BAM/CRAM file name without extension), in DIR/<name>/:
#   <name>_R1.fastq.gz, <name>_R2.fastq.gz   paired reads
#   <name>_singletons.fastq.gz               reads whose mate is missing (only written if any)
#   <name>.hg38.bam (+ .bai)                 bwa mem alignment to hg38, coordinate sorted
#   <name>.hg38.flagstat.txt, <name>.hg38.contigs.txt (reads per hg38 contig)
#   <name>.regions_hg19.bed                  regions the reads were taken from
set -euo pipefail

VARS=(); INS=(); OUT=./realigned_hg38; PAD=500; T=32; FAST=0; KEEP=0
REF19=/root/naser/liftover/ref/hg19.fa.gz
IDX38=/root/naser/liftover/ref/hg38_noalt_bwa/hg38_noalt.fa
export PATH=/data1/HUVYA/opt/miniforge3/envs/liftover/bin:/usr/local/bin:$PATH

usage() { sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }
while getopts "v:i:o:p:t:r:x:fkh" o; do
  case $o in
    v) VARS+=("$OPTARG");; i) INS+=("$OPTARG");; o) OUT=$OPTARG;; p) PAD=$OPTARG;; t) T=$OPTARG;;
    r) REF19=$OPTARG;; x) IDX38=$OPTARG;; f) FAST=1;; k) KEEP=1;; h) usage 0;; *) usage;;
  esac
done
[[ ${#VARS[@]} -gt 0 && ${#INS[@]} -gt 0 ]] || usage
for x in samtools bwa; do command -v $x >/dev/null || { echo "missing: $x" >&2; exit 1; }; done
[[ -s $IDX38.bwt ]] || { echo "hg38 bwa index not found: $IDX38.bwt" >&2; exit 1; }
[[ -s $REF19.fai ]] || { echo "hg19 reference .fai not found: $REF19.fai" >&2; exit 1; }

# -i may point to a list file
inputs=()
for i in "${INS[@]}"; do
  if [[ $i =~ \.(bam|cram)$ ]]; then inputs+=("$i"); else mapfile -t -O ${#inputs[@]} inputs < <(grep -v '^\s*$' "$i"); fi
done
mkdir -p "$OUT"
log() { echo "[$(date +%H:%M:%S)] $*"; }

# ---- hg19 regions (chr, start0, end) from all variant files -------------------------------
REG_ALL=$OUT/regions_hg19.unmerged.bed
: > "$REG_ALL"
for v in "${VARS[@]}"; do
  awk -F'\t' -v pad="$PAD" -v OFS='\t' '
    NR==1 { h=$0; sub(/^#/,"",h); n=split(h,H,"\t"); for(k=1;k<=n;k++) ix[H[k]]=k
            if ("hg19_Chr" in ix) {c=ix["hg19_Chr"]; p=ix["hg19_Pos"]; r=ix["hg19_Ref"]}
            else if ("Chr" in ix) {c=ix["Chr"]; p=ix["Pos"]; r=ix["Ref"]}
            else {c=1; p=2; r=3; if ($2 ~ /^[0-9]+$/) hdr=0; else next}
            if (("hg19_Chr" in ix) || ("Chr" in ix)) next }
    /^#/ || $p !~ /^[0-9]+$/ { next }
    { len = ($r=="_" || $r=="-") ? 1 : length($r); s=$p-1-pad; if (s<0) s=0; print $c, s, $p-1+len+pad }' "$v" >> "$REG_ALL"
done
[[ -s $REG_ALL ]] || { echo "no variants read from ${VARS[*]}" >&2; exit 1; }
log "variants: $(wc -l < "$REG_ALL") (padding ${PAD} bp)"

merge_bed() { sort -k1,1 -k2,2n | awk -v OFS='\t' '
  $1!=c || $2>e { if (c!="") print c,s,e; c=$1; s=$2; e=$3; next } { if ($3>e) e=$3 } END { if (c!="") print c,s,e }'; }

# ---- per sample --------------------------------------------------------------------------
for in in "${inputs[@]}"; do
  [[ -s $in ]] || { echo "not found: $in" >&2; continue; }
  name=$(basename "$in"); name=${name%.cram}; name=${name%.bam}
  d=$OUT/$name; mkdir -p "$d"; tmp=$d/tmp; mkdir -p "$tmp"
  log "== $name"

  # match chromosome naming (chr1 vs 1) to this file
  if samtools view -H "$in" | grep -q $'^@SQ\tSN:chr'; then addchr='chr'; else addchr=''; fi
  sed -E 's/^chr//' "$REG_ALL" | awk -v a="$addchr" -v OFS='\t' '{ $1=a $1; if ($1=="chrMT") $1="chrM"; print }' \
    | merge_bed > "$d/$name.regions_hg19.bed"
  log "regions: $(wc -l < "$d/$name.regions_hg19.bed") merged intervals"

  # read group for the realignment
  sm=$(samtools view -H "$in" | awk -F'\t' '/^@RG/{for(i=2;i<=NF;i++) if($i~/^SM:/){print substr($i,4); exit}}')
  rg="@RG\tID:${name}_hg38\tSM:${sm:-$name}\tPL:ILLUMINA\tLB:${sm:-$name}"

  # 1. names of every read (incl. secondary/MAPQ 0) touching the regions
  samtools view -@ "$T" -T "$REF19" -M -L "$d/$name.regions_hg19.bed" "$in" | cut -f1 | sort -u -S 2G -T "$tmp" > "$tmp/names.txt"
  log "read names in regions: $(wc -l < "$tmp/names.txt")"

  # 2. all primary records of those reads (mates included unless -f), 3. collate -> paired FASTQ
  if (( FAST )); then
    samtools view -@ "$T" -T "$REF19" -u -F 0x900 -M -L "$d/$name.regions_hg19.bed" "$in"
  else
    samtools view -@ "$T" -T "$REF19" -u -F 0x900 -N "$tmp/names.txt" "$in"
  fi | samtools collate -@ "$T" -u -O -T "$tmp/collate" - \
     | samtools fastq -@ "$T" -n -c 6 \
         -1 "$d/${name}_R1.fastq.gz" -2 "$d/${name}_R2.fastq.gz" \
         -s "$d/${name}_singletons.fastq.gz" -0 "$d/${name}_other.fastq.gz" - 2> "$tmp/fastq.log"
  grep -h 'processed' "$tmp/fastq.log" | sed 's/^/    /' || true
  pairs=$(( $(zcat "$d/${name}_R1.fastq.gz" | wc -l) / 4 ))
  log "read pairs: $pairs"

  # 4. realign to hg38
  if (( pairs > 0 )); then
    bwa mem -t "$T" -M -R "$rg" "$IDX38" "$d/${name}_R1.fastq.gz" "$d/${name}_R2.fastq.gz" 2> "$tmp/bwa.log" \
      | samtools sort -@ "$T" -m 1G -T "$tmp/sort" -o "$d/$name.hg38.bam" -
    samtools index -@ "$T" "$d/$name.hg38.bam"
    samtools flagstat -@ "$T" "$d/$name.hg38.bam" > "$d/$name.hg38.flagstat.txt"
    samtools idxstats "$d/$name.hg38.bam" | awk '$3>0' | sort -k3,3nr > "$d/$name.hg38.contigs.txt"
    log "hg38 BAM: $d/$name.hg38.bam  ($(awk 'NR==1{print $1}' "$d/$name.hg38.flagstat.txt") records, $(grep -m1 'primary mapped' "$d/$name.hg38.flagstat.txt" | grep -o '([^:]*' | tr -d '('))"
  fi
  for f in "$d/${name}_other.fastq.gz" "$d/${name}_singletons.fastq.gz"; do   # drop empty FASTQs
    [[ -n $(zcat "$f" | head -c1) ]] || rm -f "$f"
  done
  (( KEEP )) || rm -rf "$tmp"
done
(( KEEP )) || rm -f "$REG_ALL"
log "done -> $OUT"
