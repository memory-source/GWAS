verdict: yes
confidence: medium
integrity_status: pass (EXPERIMENT_AUDIT_DEEP.json; round-2 re-audit codex thread 01a14455-c02b-7c93-b103-24ed8aace24a reached FAIL → WARN → PASS within the declared scope)
round: deep-search direction (selection-input spectral tomography; s61/s62/s63/s64/s65/s66)

Judge: codex / gpt-6-astra (read-only, cross-family vs deepseek-flash executor), thread
01a14176-aaa4-7bf0-bf4c-9065eda1a1fe, ten rounds: score 4 → 5 → 5 → 6 → 6 → 6.5 →
6.5 → 6.5 → 7 → 7; verdict not-ready → almost → **ready** at round 10. Final:
ready (bounded measurement paper under the declared information conditions).
Evidence: results/raw/s61_spectral_tomography_chr22.json,
s61_spectral_tomography_chr1.json, s62_sample_size_chr22.json,
s63_plink_export_chr22.json, s64_saturation_chr22.json,
s65_identity_ablation_chr22.json, s66_control_provenance_chr22.json,
code/s61_spectral_tomography.py,
code/s62_sample_size.py, code/s63_plink_export.py, code/s64_saturation.py,
code/s65_identity_ablation.py, code/s66_control_provenance.py,
refine-logs/EXPERIMENT_PLAN_DEEP_20261016.md,
refine-logs/EXPERIMENT_RESULTS_DEEP_20261016.md,
review-stage/DEEP_REVIEW_20261016.md, EXPERIMENT_AUDIT_DEEP.md.

## Claim map (deep round; separate from the passive paper's CLAIMS_FROM_RESULTS.md)

| Claim | Verdict |
|-------|---------|
| D1 — the final-output identity $Y_t = V_t\Lambda_t - K_t V_t = C_H V_t$ holds on the real marker layout under the released preprocessing | yes (residual 3.8e-11 float64; 6.2e-6 through a real PLINK 2.0 export read-back) |
| D2 — an AF-preserving block swap leaves per-SNP histograms, AF and variance exactly unchanged while individual heterozygosity moves | yes (s61 witness; histograms identical, AF/var diff ≤ 4.4e-16, het shift ≤ 4.9e-3) |
| D3 — visible-subspace rank/share and scatter recovery grow with the number of final outputs; the passive strategy stays flat | yes (s61/s62/s64; monotone visible rank/share; decreasing span-restricted scatter error; passive flat) |
| D4 — accumulating several final outputs yields a membership increment over the strongest passive baseline on equal information | yes at ideal precision (q8 +0.015 EUR / +0.123 3pop / +0.013 EAS, 5/5, CI excludes zero; +0.014 through a real PLINK export, 5/5); scoped below |
| D5 — the partner's swaps leave the honest analysis's summary statistics (inflation factor, causal min-p under a simulated phenotype) essentially unchanged | yes (sci_cost; flat λ_GC and causal min-p across releases) |
| D6 — adaptive block selection beats random block swaps | no (adaptive-minus-random max +0.014, ~+0.010 or less elsewhere; CI excludes zero at chr1 3pop q4 +0.0097, chr1 EAS q4 +0.0053, chr22 EAS q4 +0.0031) |
| D7 — the scatter reconstruction survives realistic export precision | no (float64-only witness; error 1.7–4.7 at float32/4 dp) |
| D8 — the increment is specific to the tomography identity: a raw-eigenvector projection attack does not reproduce it (ablation control) | yes (s65: identity-stack q8 gap +0.015/+0.123/+0.013 vs raw-eigenvector stack −0.004/−0.012/−0.050 on the same splits; raw-V clears the baseline only at isolated EUR/3pop splits) |
| D9 — the no-release control gives no separation under three preprocessing choices (ablation control) | yes (s66: for these tested projection attacks, all 27 five-split mean AUROCs lie between 0.455 and 0.509) |

## What the results support
- D1 on real chr22/chr1 genomes (s61) and through a real PLINK 2.0 `--pca allele-wts`
  export read back from `.eigenvec.allele`/`.eigenval` (s63).
- D2: marginal QC is passed (histograms/AF/variance exact); individual statistics
  are not preserved.
- D3: coverage and reconstruction improve with q = 1/2/4/8; passive flat.
- D4: the accumulated q8 attack beats the stronger of composition-balanced
  conditioned AF and the conditioned single release; the effect appears from q≈4.
  It replicates on chr1 and grows at a larger immediate honest cohort (n=180, s62,
  one nested step, not a scaling law), and it survives the real PLINK export.
  Practical reading under the ceiling-censored AUC: TPR@5%FPR = 1.000 (EUR/EAS) and
  0.580 (3pop) at q8, against 0.916 / 0.299 / 0.899 for the strongest single-release
  baseline. The q-saturation cell (s64) shows the within-population membership
  increment saturates by q≈8 while visible share keeps rising to 0.92–0.93; the
  cross-population cohort with headroom keeps gaining (+0.049 random, +0.048
  adaptive from q8 to q16).
- D5: the honest study's joint-PCA covariate-adjusted association is essentially
  unchanged across releases.
- D6's stronger-threat reading: **uniformly random AF-preserving swaps suffice**;
  a defender cannot rely on detecting targeted input choice, which is why the
  adaptive attempt (D6) adds nothing.

## Interface scope
The claim covers interfaces that return the joint top-k loadings and eigenvalues per
re-fit (e.g. a PLINK 2.0 `--pca allele-wts` export, or SF-PCA-style selective
decryption of a study-fitted loading matrix), not all federated GWAS pipelines.

## What the results do not support
- **One fixed release never beats the strong passive baselines** (q1 gap
  −0.076/−0.129/−0.072); the whole effect is the accumulated part, not universal
  dominance.
- **The scatter recovery is a float64-only witness**; at float32/4 dp the recovered
  span gains spurious directions (relative error 1.7–4.7). Only D4's membership
  signal is claimed at export precision.
- **No adaptive advantage** (D6).
- **Small scope**: real 1000G chr22/chr1 prefixes; honest n=120/180/360; one k=10,
  block length 20, candidate budget 24; q ≤ 8 in the headline cells (q = 16 only in
  the within:EUR/3pop/EAS chr22 saturation cell); one software package (PLINK 2.0),
  not a full federated pipeline; roles drawn from the same 1000G release.
- No honest-cohort reconstruction, no family identification, no real phenotype;
  the GWAS rows use a simulated phenotype.
