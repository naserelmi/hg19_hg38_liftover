# Iranome hg19 → hg38 liftover check

Run on NPU: `/root/naser/liftover` (scripts: `run.sh`, `annovar2vcf.py`, `liftvcf.py`, `compare.py`, `provided_raw.py`; outputs in `work/`).
References: UCSC hg19.fa.gz and hg38.fa.gz, UCSC `hg19ToHg38` / `hg38ToHg19` over.chain, UCSC `liftOver -minMatch=1`.

## Method
1. hg19 table (ANNOVAR style) → VCF with anchor bases from hg19. All 1,594,011 Ref alleles match hg19.
2. `bcftools norm -f hg19` (left-align and normalize): 77,986 indels were realigned.
3. Lift hg19 → hg38: lift the REF span with liftOver, re-read REF from hg38, reverse-complement on `-` chains, then `bcftools norm -f hg38`.
4. Lift back hg38 → hg19 the same way, normalize, and compare with step 2 per variant (ID = original row).
5. Audit the provided `hg38_iranome.txt` against the correct liftover (row level, matched by lifted Start and AF).

CrossMap 0.7.4 was tried first and dropped. It rejects ~90% of deletions as "REF==ALT" and truncates the REF of lifted deletions.

## A. Round trip hg19 → hg38 → hg19 (correct liftover)

| status | n | % |
|---|---:|---:|
| IDENTICAL (same chr, pos, ref, alt) | 1,585,113 | 99.44 |
| failed hg19→hg38 | 5,307 | 0.33 |
| changed position and allele (other paralog, opposite strand) | 1,597 | 0.10 |
| changed position (same chr) | 1,341 | 0.08 |
| changed chromosome | 355 | 0.02 |
| failed hg38→hg19 | 298 | 0.02 |

Reasons for failing hg19→hg38: 2,134 REF/ALT swapped (hg38 reference = hg19 ALT, so the AF would need inverting), 1,287 deleted in hg38, 573 partially deleted, 1,172 landing only on alt/random contigs, 127 span length changed, 14 split.
Another 171 variants lift fine, but the hg38 REF base differs from hg19 (`REFDIFF` tag).
Changed variants cluster in segmental duplications (chr1 1q21/1p36, chr9 pericentromeric, chr22q11, chr10). Example: chr1:1582880 G>C maps to hg38 chr1:1714674, and that maps back to chr1:1646113 (the other CDK11 copy).

Files: `work/rt.roundtrip_all.tsv.gz` (every variant), `work/rt.roundtrip_changed.tsv` (the 8,898 non-identical ones), `work/fwd38.unmap`, `work/back19.unmap`.

## B. Is the provided hg38_iranome.txt correct? — No

Of the 1,594,011 hg19 rows:

| category | n |
|---|---:|
| identical to a correct liftover | 1,049,996 (65.9%) |
| SNV → 1-bp deletion (ALT allele lost, e.g. `1:13079 C>G` → `C>-`) | 389,100 |
| deletion truncated to its first base (e.g. `GACA>-` → `G>-`) | 38,150 |
| insertion turned into a deletion of the inserted sequence | 25,350 |
| insertion anchored wrongly: anchor base overwrote the first inserted base (`->GAGA` → `G>GAGA`, should be `G>GGAGA`) | 22,952 |
| insertion → SNV / other garbled alleles | 5,884 |
| minus-strand regions: REF reverse-complemented but ALT not (`C>T` → `G>T`, should be `G>A`) | ~3,900 |
| missing from the provided file | 56,109 |
| not liftable | 1,859 |

The provided file also has 20,917 malformed rows where the Alt column holds AF values (`A  0.0089,0.62  .`). It has 27,511 rows with no hg19 counterpart.
**Conclusion:** coordinates were lifted correctly, but alleles were corrupted for about one third of the variants. The provided hg38 file should not be used. Regenerate it from hg19 (`work/fwd38.norm.vcf`).

## C. Issues in the original hg19 file
* 10,767 exact duplicate keys (same chr/start/ref/alt, different AF), e.g. `17:7349623 ->TGTGTG` AF 0.0113 and 0.0127.
* After normalization there are 25,413 duplicate keys (equivalent indel representations in repeats, e.g. `ACACACACAC>AC` = `ACACACAC>-`).
