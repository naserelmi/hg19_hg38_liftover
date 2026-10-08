# DataIranome hg19 → hg38 liftover

Lifts the Iranome allele-frequency table from hg19 to hg38, checks every variant with an hg38 → hg19 round trip, and sorts the variants that don't come back identical by reason.
Scripts use absolute paths for the NPU set-up (`/root/naser/liftover`, conda env `liftover`); change `D=` / `PATH` in `code/run_pipeline.sh` and the defaults in `code/realign_unmatched.sh` for another machine. Data, references and results are not in this repository.

### Set-up
    # references (UCSC)
    mkdir -p ref && cd ref
    for u in hg19/bigZips/hg19.fa.gz hg38/bigZips/hg38.fa.gz hg19/liftOver/hg19ToHg38.over.chain.gz hg38/liftOver/hg38ToHg19.over.chain.gz; do
      wget https://hgdownload.soe.ucsc.edu/goldenPath/$u; done
    for g in hg19 hg38; do zcat $g.fa.gz | bgzip > $g.bgz.fa.gz && mv $g.bgz.fa.gz $g.fa.gz && samtools faidx $g.fa.gz; done
    # hg38 no-alt bwa index (for realign_unmatched.sh)
    mkdir -p hg38_noalt_bwa && cut -f1 hg38.fa.gz.fai | grep -vE "_alt|_fix" > hg38_noalt_bwa/contigs.txt
    samtools faidx hg38.fa.gz -r hg38_noalt_bwa/contigs.txt > hg38_noalt_bwa/hg38_noalt.fa && bwa index hg38_noalt_bwa/hg38_noalt.fa
    # tools
    mamba create -n liftover -c conda-forge -c bioconda python=3.11 pysam ucsc-liftover   # + bcftools, samtools, bwa

### Run
From `/root/naser/liftover`:

    bash code/run_pipeline.sh > pipeline.log 2>&1

Takes about 5 minutes. Intermediate files go to `work_data/` and result tables to `results/`.

## Inputs
* `DataIranome.txt`: hg19, no header. Columns are chr, pos, ref, alt, het, hom, AF. `_` marks an empty allele, using the ANNOVAR indel convention.
  AF = (het + 2·hom) / AN. AN is recovered per site and is at most 1600 (800 people).
* `ref/`: UCSC hg19.fa.gz and hg38.fa.gz (bgzip + faidx), plus the UCSC hg19ToHg38 and hg38ToHg19 over.chain.gz files.
* conda env `liftover`, which provides python3, pysam and ucsc-liftover; bcftools comes from /usr/local/bin.

## Code (`code/`)
| script | role |
|---|---|
| `run_pipeline.sh` | runs the whole pipeline |
| `annovar2vcf.py` | builds a VCF from the table: adds anchor bases from the genome, checks Ref against the genome, and keeps het/hom/AC/AN in INFO |
| `liftvcf.py` | lifts the REF span with `liftOver -minMatch=1`, re-reads REF from the target genome, reverse-complements on `-` chains, and flags ref/alt swaps, alt contigs and span changes |
| `compare.py` | classifies each variant after the hg19 → hg38 → hg19 round trip |
| `make_outputs.py` | writes the six result tables and `summary.txt` |
| `vcf2table.py` | converts a normalized VCF back to a DataIranome-style table |

Every VCF is sorted, left-aligned and normalized with `bcftools norm -f <genome>`.
CrossMap 0.7.4 was not used, because its `vcf` mode drops and truncates deletions.

## Results (`results/`)
Rows are in DataIranome style: `_` = empty allele, a deletion's Pos is its first deleted base, and an insertion's Pos is the base before the insertion.

| file | n | meaning |
|---|---:|---|
| 1_matched_variants_hg38.tsv | 1,604,163 | lifting back gives the identical hg19 variant. `hg38_ref_differs_from_hg19`=yes for 107 |
| 2_altswapped_variants_hg38.tsv | 1,734 | the hg38 reference base is the hg19 ALT, so REF/ALT are swapped. Biallelic sites (1,420): Het unchanged, Hom = N − Het − Hom, AF = 1 − AF. Multiallelic sites (314): AC = AN − AC − AC(other alts), and genotypes are NA. 521 end up with AC = 0 (all 800 were hom-alt in hg19), so they are not variants in hg38 |
| 3_deleted_variants_hg38.tsv | 1,134 | the region is deleted (801) or partially deleted (333) in hg38, so there is no hg38 position |
| 4_chrom_changed_variants_hg38.tsv | 290 | lifts to hg38, but lifting back lands on another chromosome (segmental duplication) |
| 5_broken_variants_hg38.tsv | 663 | maps only to alt/unplaced contigs (434), span length changed (69), split (11), or the hg38 position does not lift back (149) |
| 6_duplicate_regions_mismatch_hg38.tsv | 579 | lifts to hg38, but lifting back lands elsewhere on the same chromosome (`shift_bp`) |

Other files:
* `all_lifted_variants_hg38.tsv` and `.vcf.gz` contain every variant with an hg38 position. These are files 1 + 4 + 6 plus the 149 not-liftable-back variants from file 5.
* `roundtrip_per_variant.tsv.gz` gives each variant's hg19, hg38 and round-trip result.

## Realigning reads behind unmatched variants: `code/realign_unmatched.sh`
`results/unmatched_variants_hg19.tsv` holds all 4,400 variants from tables 2–6, with category and reason.
For each hg19 BAM or CRAM, the script:
1. takes every read touching ±500 bp around those variants, including secondary and MAPQ 0 reads;
2. scans the whole file for those read names, so mates outside the regions are included too;
3. writes `<name>_R1.fastq.gz` and `<name>_R2.fastq.gz`;
4. realigns them with `bwa mem -M` to hg38 no-alt (`ref/hg38_noalt_bwa/`), giving `<name>.hg38.bam`.

    code/realign_unmatched.sh -v results/unmatched_variants_hg19.tsv -i sample.cram -o realigned_hg38 -t 32
    code/realign_unmatched.sh -v results/6_duplicate_regions_mismatch_hg38.tsv -i cram_list.txt   # one category, many samples

`-p` sets the padding, `-f` skips the whole-file mate scan (faster), `-x` sets another hg38 bwa index, and `-r` sets the hg19 reference used to decode CRAMs. `-h` prints help.
On one WES CRAM (RBG1037, about 0.9 GB) it took about 2 minutes: 321k read pairs, no orphan mates, 99.9% mapped to hg38.
