# External data

The genotype VCFs are **not** tracked in this repository (2.5 GB). Everything below can
be fetched from public sources.

## 1000 Genomes Phase 3

Expected filenames under `data/external/`:

| file | content |
|---|---|
| `chr22.vcf.gz` | chromosome 22 genotypes |
| `ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | chromosome 1 |
| `ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | chromosome 2 |
| `1000g_panel.txt` | sample → population / super-population table (tracked) |

Source: the 1000 Genomes Project Phase 3 release, `phase3_shapeit2_mvncall_integrated_v5b.20130502`
(`ftp.1000genomes.ebi.ac.uk`, path `/vol1/ftp/release/20130502/`). The chromosome 22 file
is the same release, renamed to `chr22.vcf.gz`.

The panel table has 2504 samples with the columns
`sample, pop, super_pop, gender`; the super-population counts are
AFR 661, AMR 347, EAS 504, EUR 503, SAS 489.

Loaders: `code/real_panel.py` and `code/s48_real_structure.py` read these paths
(`ROOT/data/external/...`). Parsed panels are cached under `results/cache/` (not tracked).

## Genetic maps

`data/external/genetic_map/genetic_map_chr{1,2,22}.txt` are the 1000G Phase-3 IMPUTE2
genetic maps: whitespace-delimited `position  COMBINED_rate(cM/Mb)  Genetic_Map(cM)`.
They are small (~22 MB total) and are not tracked; fetch them from the IMPUTE2 map
distribution referenced by the 1000G release. `code/genetic_map.py` parses them with
`numpy.loadtxt` and linearly interpolates cumulative cM at an arbitrary bp.

## PLINK 2.0

`code/s63_plink_export.py` and `code/s72_printed_matched.py` invoke a `plink2` binary at
`scratch/plink2bin/plink2` (not tracked). Any PLINK 2.0 build works (the results quoted in
the paper were produced with `PLINK v2.0.0-a.7.13LM AVX2 Intel`). The scripts write a
minimal VCF, run `--pca <k> allele-wts [approx]`, and read the printed
`.eigenvec.allele` (ALT rows) and `.eigenval` back.
