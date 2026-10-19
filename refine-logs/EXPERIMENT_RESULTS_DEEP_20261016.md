# Experiment Results — Selection-Input Spectral Tomography (deep-search round, 2026-10-16)

Script: `code/s61_spectral_tomography.py` → `results/raw/s61_spectral_tomography_chr22.json`
(one run, 135–547 s depending on machine load, 5 split seeds per cohort, chromosome 22,
1.5 kb thinning). (Second chromosome run: `results/raw/s61_spectral_tomography_chr1.json`.)
Ablation controls: `code/s65_identity_ablation.py` → `results/raw/s65_identity_ablation_chr22.json`;
`code/s66_control_provenance.py` → `results/raw/s66_control_provenance_chr22.json`.
The q-saturation cell is `code/s64_saturation.py` → `results/raw/s64_saturation_chr22.json`.
Honest cohort = real 1000G chr22 genotypes. Marker counts: `within:EUR` 3,743;
`super:3pop` 3,271; `within:EAS` 3,271. The loader takes the first 20,000
MAF-filtered variants before thinning, so these are chromosome-22 prefixes
(≈16.05–24.56 / 16.05–23.22 Mb), not chromosome-wide samples.

Roles per split (disjoint, matched per group): honest / partner / non-member /
attacker reference / shadow with per-group caps 120 / 120 / 80 / 60 / 40. Honest n
is the pooled count (120 for the within-population arms, 360 for `super:3pop`).
Roles are drawn without any phenotype; the phenotype is simulated and used only by
the GWAS-power row, whose association is fitted on the honest study alone.

Output accounting: output 1 is the shared passive baseline release for every
strategy; each further output always starts from a fresh partner input, so output q
costs q releases for every strategy.

## B0 — Interface witness + preservation audit (real marker layout, real LD)
- `max_eigen_equation_residual` = 3.8e-11.
- `histograms_identical_all_releases` = true everywhere. Exact per-SNP histogram
  equality over all 12 releases implies the observed-genotype-count-based HWE
  statistic and the per-SNP missing count are also unchanged (both are functions of
  the per-SNP genotype counts alone); the block swap moves no allele between rows.
- `max_abs_af_change` = 0.0; `max_abs_var_change` = 4.4e-16.
- `genotypes_changed_per_release` (matrix entries changed; two per differing
  marker) ranges 0 (a matched block is a no-op) to 40.
- `max_abs_row_heterozygosity_change` = 1.8e-3 – 4.9e-3 (the individual-level
  statistic that does move).

## Baselines under the SAME allowed information
The attacker is granted the released AF, the partner matrix and the independent
shadow panel (neither contains the scored targets). Held-out-half oriented:
`af_ref` (background = reference alone), `af_pooled` (background pooled over
reference+partner+shadow), `af_cond` (composition-balanced removal of the known
calibration members from the released AF, pooled background), `spec_cond` (the
single passive release with the calibration members' scatter contribution removed,
scored by its range; orientation fixed on the unconditioned release).

| cohort | af_ref | af_pooled | af_cond | spec_cond |
|---|---:|---:|---:|---:|
| within:EUR | 0.689 | 0.791 | 0.882 | 0.985 |
| super:3pop | 0.633 | 0.695 | 0.784 | 0.757 |
| within:EAS | 0.688 | 0.812 | 0.915 | 0.987 |

`gapStrong` = paired increment over the stronger of `af_cond` and `spec_cond`,
with a 10k-resample bootstrap 95% CI. `TPR@5%FPR` is the held-out-half
true-positive rate at a 5% false-positive rate (`common.tpr_at_fpr`); because the
stacked attack reaches AUC 1.000, the AUC gap understates the practical effect.

## No-release control (the "isn't this just knowing some members?" objection)
The attacker is granted the KNOWN calibration members' full genotypes but produces
NO release, and scores membership by projection onto (i) their own top-k PCA span,
(ii) their top-(qk) span, (iii) a random k-dim span. Held-out oriented. The rows
entering the control score are bounded at |z| ≤ 5: the `super:3pop` and `within:EAS`
panels are MAF-filtered on the full 2504-sample 1000G panel, so they carry
globally-common markers that are monomorphic inside the cohort (cohort `mu == 0`,
392 of 3,271 chr22-EAS markers); without the bound the scale clip (1e-6) gives an
out-of-`mu` carrier |z| up to 2×10³ and a random span picks up a spurious
anti-membership signal (the unclipped EAS `random_k` read 0.579, an artifact — see
the review log). When the |z| ≤ 5 bound is applied it changes **all three** control
arms, since the bound is applied to the rows entering the control score; the effect
is largest on the random-span arm (EAS `random_k` 0.579 → 0.463) and small on the
cal-PCA arms (EAS `cal_pca_k` 0.490 → 0.497, `cal_pca_qk` 0.450 → 0.455), which
already have near-zero weight at `mu == 0` markers. The attack scores themselves are
not clipped (the attack uses the released operator, not a random span), and every
headline number is unchanged.

| cohort | cal PCA k=10 | cal PCA qk=80 | random k=10 |
|---|---:|---:|---:|
| within:EUR | 0.467 | 0.459 | 0.488 |
| super:3pop | 0.501 | 0.494 | 0.485 |
| within:EAS | 0.497 | 0.455 | 0.463 |

All at or just above chance (largest 3pop `cal_pca_k` 0.501): the release-derived span
carries the membership signal, not the calibration-member genotypes. The single-release
TPR@5%FPR is 0.675 (EUR), 0.139 (3pop), 0.616 (EAS); the conditioned single release
reaches 0.916 / 0.299 / 0.899.

## B1/B2/B3 — Coverage and membership vs number of final outputs

### within:EUR
| strategy | q | rank_vis | visible_share | scatter_recon_err | AUC | TPR@5%FPR | gapStrong |
|---|---:|---:|---:|---:|---:|---:|---:|
| passive  | 1 | 10 | 0.593 | 0.805 | 0.908 | 0.675 | −0.076 [−0.094,−0.058] |
| passive  | 8 | 10 | 0.593 | 0.805 | 0.908 | 0.675 | −0.076 [−0.094,−0.058] |
| random   | 2 | 18 | 0.661 | 0.749 | 0.961 | 0.837 | −0.024 [−0.066,+0.001] |
| random   | 4 | 37 | 0.761 | 0.648 | 0.996 | 0.981 | +0.011 [+0.003,+0.020] |
| random   | 8 | 72 | 0.853 | 0.521 | 1.000 | 1.000 | +0.015 [+0.009,+0.022] |
| adaptive | 8 | 75 | 0.860 | 0.510 | 1.000 | 1.000 | +0.015 [+0.009,+0.022] |

### super:3pop
| strategy | q | rank_vis | visible_share | scatter_recon_err | AUC | TPR@5%FPR | gapStrong |
|---|---:|---:|---:|---:|---:|---:|---:|
| passive  | 1 | 10 | 0.920 | 0.393 | 0.655 | 0.139 | −0.129 [−0.141,−0.115] |
| passive  | 8 | 10 | 0.920 | 0.393 | 0.655 | 0.139 | −0.129 [−0.141,−0.115] |
| random   | 2 | 20 | 0.930 | 0.367 | 0.739 | 0.213 | −0.045 [−0.060,−0.028] |
| random   | 4 | 39 | 0.943 | 0.333 | 0.830 | 0.384 | +0.047 [+0.031,+0.063] |
| random   | 8 | 75 | 0.956 | 0.294 | 0.907 | 0.580 | +0.123 [+0.103,+0.141] |

### within:EAS
| strategy | q | rank_vis | visible_share | scatter_recon_err | AUC | TPR@5%FPR | gapStrong |
|---|---:|---:|---:|---:|---:|---:|---:|
| passive  | 1 | 10 | 0.602 | 0.798 | 0.915 | 0.616 | −0.072 [−0.082,−0.062] |
| passive  | 8 | 10 | 0.602 | 0.798 | 0.915 | 0.616 | −0.072 [−0.082,−0.062] |
| random   | 2 | 18 | 0.667 | 0.743 | 0.964 | 0.810 | −0.023 [−0.045,−0.010] |
| random   | 4 | 37 | 0.768 | 0.640 | 0.991 | 0.956 | +0.004 [+0.001,+0.006] |
| random   | 8 | 70 | 0.854 | 0.520 | 0.999 | 1.000 | +0.013 [+0.010,+0.015] |
| adaptive | 8 | 75 | 0.866 | 0.501 | 0.999 | 1.000 | +0.012 [+0.010,+0.015] |

## B5 — Export-precision probe (random q8 stack, released V/Lambda re-formed at each
precision; `gapStrong` uses a PRECISION-MATCHED strong passive baseline)
| cohort | precision | scatter_recon_err | AUC | gapStrong |
|---|---|---:|---:|---:|
| within:EUR | float64 | 0.521 | 1.000 | +0.0151 |
| within:EUR | float32 | 1.811 | 1.000 | +0.0151 |
| within:EUR | 6 dp | 2.099 | 0.996 | +0.0115 |
| within:EUR | 4 dp | 2.620 | 0.963 | −0.0213 |
| super:3pop | float64 | 0.294 | 0.907 | +0.1233 |
| super:3pop | float32 | 2.895 | 0.893 | +0.1096 |
| super:3pop | 6 dp | 4.123 | 0.806 | +0.0224 |
| super:3pop | 4 dp | 4.694 | 0.678 | −0.1061 |
| within:EAS | float64 | 0.520 | 0.999 | +0.0126 |
| within:EAS | float32 | 1.708 | 0.997 | +0.0106 |
| within:EAS | 6 dp | 2.172 | 0.990 | +0.0028 |
| within:EAS | 4 dp | 2.683 | 0.953 | −0.0336 |

## B4 — Scientific-cost row (simulated phenotype, 8 causal markers)
Two parts:
- Honest study alone with its own top-10 PCs: `n_study` 120 (within arms) / 360
  (3pop); `lambda_gc` 0.90 – 1.52; `min_p_causal` 5.8e-7 – 0.107; `n_causal_in_mask`
  3–8 of 8.
- Per-release joint-PCA covariates (`sci_cost`): the honest study's association is
  re-fit using the covariates from each round's JOINT PCA of (honest, partner_t).
  Over the eight random-strategy releases, the mean `lambda_gc` is 1.09–1.10 (EUR),
  1.15 (3pop), 1.24–1.26 (EAS) and the mean causal min-p is 0.026 (EUR), 0.020
  (3pop), 0.035 (EAS), essentially flat across releases: the partner's AF-preserving
  swaps do not materially perturb the honest study's covariate-adjusted analysis.

## What these results support
1. The identity $Y_t = C_H V_t$ holds to $3.8\times10^{-11}$ on the real marker
   layout (B0).
2. An AF-preserving block swap leaves per-SNP histograms, AF and variance exactly
   unchanged (individual heterozygosity shifts up to 4.9e-3): marginal QC is
   passed, individual structure is not preserved (B0).
3. Visible rank and visible share increase monotonically with the number of final
   outputs q=1/2/4/8; the scatter reconstruction error decreases monotonically; the
   passive strategy stays flat throughout. (B1–B3)
4. **Under ideal float64 precision, accumulating several final outputs beats the
   strongest passive baseline on equal information** (the stronger of conditioned
   AF and the conditioned single release): +0.015 EUR, +0.123 3pop, +0.013 EAS at
   q8, 5/5 splits, CI excluding zero, on both chromosomes. At q4 the accumulated
   increment's CI excludes zero on chr22 for all three cohorts (5/5 positive splits
   except EUR 4/5) and on chr1 for EUR only; the within-population increments are
   small. The membership quantity is a study-inclusion score; under the
   ceiling-censored AUC at 1.000 the practical reading is TPR@5%FPR: **1.000
   (EUR/EAS, all members recovered) and 0.580 (3pop) at q8, against 0.916 / 0.299 /
   0.899 for the strongest single-release baseline.**
5. The partner's swaps leave the honest study's covariate-adjusted summary
   statistics essentially unchanged: the joint-PCA covariate inflation factor and the
   causal-marker minimum p-value are flat across releases (`sci_cost`: `lambda_gc`
   within 1.09–1.26 of the honest-only baseline, mean causal min-p 0.020–0.035), so
   the leakage accrues without a visible change in those association summaries. This
   is a summary-statistic stability statement (inflation factor and minimum p-value
   under a simulated phenotype), not a power measurement at a fixed threshold.
6. The identity and the accumulated-output increment reproduce through a REAL PCA
   export: PLINK 2.0 writes the joint allele weights and eigenvalues, the attacker
   reads them back, and the q8 stacked attack beats the strongest single release by
   +0.014 [+0.007,+0.021] (5/5) at PLINK's printed precision.
7. **The increment is not reproduced by a raw-eigenvector projection attack.** Stacking
   the raw released eigenvectors $V_t$ (without forming $Y_t=V_t\Lambda_t-K_tV_t$)
   leaves the q8 gap at or below the strongest passive baseline in all three chr22
   cohorts on the same split seeds and swap sequences (−0.004 EUR [4/5 splits, CI
   includes zero], −0.012 3pop [2/5, CI includes zero], −0.050 EAS [0/5],
   `results/raw/s65_identity_ablation_chr22.json`), against the identity stack's
   +0.015/+0.123/+0.013 (5/5 each). The raw-$V$ gap is most negative at q1
   (−0.17/−0.19/−0.19) and improves as q grows but does not clear the baseline except
   at isolated EUR/3pop splits. This shows the increment is specific to measuring the
   honest scatter through the released object, not a property of the released span
   alone; it establishes that for this alternative attack, not necessity against all
   possible attackers.
8. **The no-release control gives no separation under each of three preprocessing
   choices (a control on the control).** The attacker holds the known calibration
   members' full genotypes but no release; the calibration PCA spans are fitted and
   scored in each of three self-consistent coordinate systems (the pipeline's public
   baseline; the calibration members' own per-marker center/scale; raw centered) —
   `results/raw/s66_control_provenance_chr22.json`. For these tested projection
   attacks, all 27 five-split mean AUROCs lie between 0.455 and 0.509 (largest upward
   deviation from 0.5 is +0.009): the calibration-genotype grant does not separate
   members under any tested transformation. This is a bound over the tested transforms
   with 5 splits, not an equivalence proof, and the attack scores themselves are
   unclipped while the two standardized control spaces are |z|-bounded (see the
   no-release-control section).
9. The scope of the interface is declared: the claim covers interfaces returning the
   joint top-k loadings and eigenvalues per re-fit (PLINK 2.0 `--pca allele-wts`;
   SF-PCA-style selective decryption), not all federated GWAS pipelines.

## What these results do NOT support
- **The single passive release never beats the strong passive baselines.** At q1
  `gapStrong` is −0.076 EUR / −0.129 3pop / −0.072 EAS: one fixed release loses to
  the conditioned spectral baseline. The whole effect is the *accumulated* part,
  and it is not a universal dominance.
- **The membership increment disappears under the coarsest rounding.** With a
  precision-matched passive baseline, honest members still separate at float32
  (gap +0.015/+0.110/+0.011), the cross-population gap falls at 6 dp (+0.022) and
  the within-population gaps approach zero (EUR +0.0115, EAS +0.0028), and at 4 dp
  every arm loses to the strong passive baseline. A 4-decimal rounding (as opposed
  to PLINK's 6-significant-digit export, s63) is therefore a boundary of the effect;
  a real 6-significant-digit export does keep a small positive increment.
- **The scatter recovery is a float64-only witness.** The relative scatter
  reconstruction error — `‖Y P_Q − C_H‖ / ‖C_H‖`, i.e. recovery of the
  span-restricted $C_H P_Q$ along the released directions, never the full $C_H$ —
  is 0.29–0.53 in float64 but 1.7–4.7 at float32 or 4 dp; the recovered span gains
  spurious directions under rounding. Scatter-recovery statements are limited to
  ideal precision and to the span-restricted operator.
- **Adaptive selection gives no material advantage over random block swaps.** The
  adaptive-minus-random increment never exceeds +0.014 (chr22 EUR q2, +0.0138, CI
  includes zero); the next largest is chr22 EAS q2 +0.0102 (CI includes zero) and
  everything else is ≤ +0.010. It excludes zero at three cells only: chr1 3pop q4
  (+0.0097 [+0.003,+0.018], gone by q8), chr1 EAS q4 (+0.0053 [+0.002,+0.009]) and
  chr22 EAS q4 (+0.0031 [+0.0005,+0.005]). No cell exceeds +0.014, so it is not a
  material advantage.
- **No honest-cohort reconstruction, no family identification, no real phenotype.**
  Membership is a study-inclusion score; the GWAS-power rows use a simulated
  phenotype. Per-release covariate perturbation is measured (`sci_cost`) as the
  change in the honest study's covariate-adjusted association summaries across the
  joint-PCA releases (not a direct PC-geometry change).
- **One chromosome prefix per run.** The chr22 run is replicated on chr1 (see below)
  and the honest-sample-size sensitivity is measured at n=120 vs n=180 (s62); the
  real-software export check uses PLINK 2.0 on chr22 EUR (s63). Cells are 120/360
  real rows on a chr prefix; one k=10, block length 20, candidate budget 24.
- **Exported result files are read into a snapshot.** The export check reads the
  objects back and re-runs the attack; the export precision (6 significant digits)
  gives a maximum observed identity residual of 6.2e-6 over the 40 exports, and the
  s63 increment (≈0.014) is small because the strongest single release is already
  near 0.99. The check reuses the public center/scale and the released AF the plan
  already grants; the two PCA files are not by themselves the whole input. It uses
  one software package (PLINK 2.0 `--pca allele-wts`), not a full federated
  pipeline.

## Chromosome-1 replication (`results/raw/s61_spectral_tomography_chr1.json`)
Run on chr1 (EUR 3,934 markers; 3pop/EAS 3,230), same protocol, 5 splits. `TPR@5%FPR`
is the random q8 arm; the parenthesised value is the strongest single-release
baseline (conditioned single release).

| cohort | af_ref | af_cond | spec_cond | random q4 gapStrong | random q8 gapStrong | q8 TPR@5%FPR (baseline) |
|---|---:|---:|---:|---:|---:|---:|
| within:EUR | 0.720 | 0.906 | 0.965 | +0.027 [+0.015,+0.039] | +0.034 [+0.021,+0.046] | 1.000 (0.861) |
| super:3pop | 0.640 | 0.805 | 0.748 | +0.013 [−0.015,+0.041] | +0.089 [+0.065,+0.114] | 0.475 (0.240) |
| within:EAS | 0.729 | 0.933 | 0.980 | +0.004 [−0.001,+0.008] | +0.016 [+0.011,+0.023] | 0.997 (0.900) |

The ordering replicates on an independent chromosome: the single release loses to
the conditioned baselines and the accumulated q8 attack beats them in all three
cohorts (5/5 splits each; q4 positive on 5/5 EUR, 3/5 3pop, 4/5 EAS).

## q-saturation cell (`results/raw/s64_saturation_chr22.json`)
The capability question — how many independent outputs the partner can design — was
tested at one design point (q ≤ 8). This cell extends the SAME 5 split seeds and the
SAME machinery to q = 16 on the three chr22 cohorts. The within-population arms are
already ceiling-censored at q8 (AUC 1.000, TPR@5%FPR 1.000) so a flat q8→q16 gap
there is partly "cannot exceed perfect"; `super:3pop` sits at AUC 0.907 / TPR@5%FPR
0.580 at q8 and has headroom.

| cohort | strategy | q8 gapStrong | q16 gapStrong | q8 visible_share | q16 visible_share | q8→q16 TPR@5%FPR |
|---|---|---:|---:|---:|---:|---:|
| within:EUR | random | +0.0151 [+0.0086,+0.0217] | +0.0153 [+0.0088,+0.0219] | 0.853 | 0.927 | 1.000 → 1.000 |
| within:EUR | adaptive | +0.0152 [+0.0087,+0.0218] | +0.0153 [+0.0088,+0.0219] | 0.860 | 0.932 | 1.000 → 1.000 |
| super:3pop | random | +0.1233 [+0.1026,+0.1414] | +0.1725 [+0.1555,+0.1880] | 0.956 | 0.969 | 0.580 → 0.775 |
| super:3pop | adaptive | +0.1277 [+0.1078,+0.1446] | +0.1755 [+0.1571,+0.1917] | 0.957 | 0.970 | 0.609 → 0.801 |
| within:EAS | random | +0.0126 [+0.0099,+0.0154] | +0.0130 [+0.0104,+0.0156] | 0.854 | 0.923 | 1.000 → 1.000 |
| within:EAS | adaptive | +0.0124 [+0.0096,+0.0152] | +0.0131 [+0.0105,+0.0157] | 0.866 | 0.934 | 1.000 → 1.000 |

The within-population membership increment saturates at q ≈ 8 (the attack is already
at 1.000 TPR@5%FPR) while the visible share keeps rising to 0.92–0.93 and the
scatter reconstruction error keeps falling (0.52 → 0.38 EUR, 0.29 → 0.25 3pop). The
cohort with headroom does NOT saturate: `super:3pop` gains a further +0.049 (random)
/ +0.048 (adaptive) in the q8→q16 gap and +0.20 in TPR@5%FPR. Saturation of the
membership effect is therefore a within-population statement, not a general one.

## Honest-sample-size sensitivity (`results/raw/s62_sample_size_chr22.json`)
Nested design: partner / nonmember / reference / shadow are IDENTICAL across the two
arms (shared 120/80/60/40 roles), and the n=120 honest rows are a subset of the
n=180 honest rows; only the honest cohort grows (chr22, 5 splits). `gapStrong` =
paired increment over the stronger of conditioned AF and the conditioned single
release.

| n_honest | strategy | af_cond | spec_cond | q4 gapStrong | q8 gapStrong | q8 TPR@5%FPR |
|---:|---|---:|---:|---:|---:|---:|
| 120 | random | 0.882 | 0.985 | +0.011 [+0.003,+0.020] | +0.015 [+0.009,+0.022] | 1.000 |
| 120 | adaptive | 0.882 | 0.985 | +0.013 [+0.007,+0.020] | +0.015 [+0.009,+0.022] | 1.000 |
| 180 | random | 0.845 | 0.948 | +0.034 [+0.021,+0.050] | +0.051 [+0.034,+0.070] | 0.995 |
| 180 | adaptive | 0.845 | 0.948 | +0.038 [+0.028,+0.050] | +0.051 [+0.034,+0.070] | 0.993 |

The single release still loses at q1 in both arms (−0.076 / −0.107). The
accumulated-output increment is larger at the larger honest cohort (q8 +0.051 vs
+0.015). This is one size step on one chromosome prefix; it does not establish a
scaling law.

## Real-software export interface (`results/raw/s63_plink_export_chr22.json`)
PLINK 2.0 (`PLINK v2.0.0-a.7.13LM`, AVX2) runs the joint PCA on the real EUR chr22
panel written as VCF; the attacker reads back `.eigenvec.allele` (allele weights,
one row per allele) and `.eigenval` and forms $Y_t=V_t\Lambda_t-K_tV_t$ in PLINK's
own variance-standardized scaling ($x=(g-2p)/\sqrt{2p(1-p)}$, eigenvalues =
scatter-eigenvalues / marker count; verified that PLINK eigenvalues equal
scatter-eigenvalues divided by $M$). Marker order matches the VCF at every marker.

| quantity | value |
|---|---|
| identity residual $\max\|C_H V_t - Y_t\| / \max\|C_H V_t\|$ (max over 40 exports) | 6.2e-6 |
| marker order check | true at every marker |
| reference-AF baseline | 0.615–0.775 |
| composition-balanced conditioned AF | 0.855–0.924 |
| single release (read back) | 0.879–0.951 |
| PLINK-conditioned single release | 0.976–0.994 |
| stacked q8 attack (read back) | 0.994–1.000 |
| paired q8 increment over the strongest passive baseline | +0.014 [+0.007,+0.021], 5/5 |
| q8 TPR@5%FPR (read back) vs strongest single release | 1.000 / 0.962 vs 0.833–0.984 |

The identity and the accumulated-output increment reproduce through a real PCA
export, at the printed 6-significant-digit precision; the residual is the maximum
observed over these 40 exports at that precision, not a float64 memory (s61's
genuine memory is 3.8e-11). The effect size is small (≈0.014) because the strongest
single release is already near 0.99 here. The read-back attack reuses the public
center/scale and the released AF that the plan already grants; the two PCA files
are not by themselves the whole input.
