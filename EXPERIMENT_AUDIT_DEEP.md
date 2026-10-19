# Experiment Audit Report — Selection-Input Spectral Tomography (deep-search round)

**Date**: 2026-10-16
**Auditor**: Codex / gpt-6-astra, read-only (cross-family; thread `01a14126-7f48-7912-9d7f-88dc574296b9`; round-2 re-audit thread `01a14455-c02b-7c93-b103-24ed8aace24a`)
**Project**: /data/lvvyongji1/GWAS
**Artifacts audited**: `code/s61_spectral_tomography.py`, `code/s62_sample_size.py`,
`code/s63_plink_export.py`, `code/s64_saturation.py`, `code/s65_identity_ablation.py`,
`code/s66_control_provenance.py`,
`results/raw/s61_spectral_tomography_chr22.json` / `_chr1.json`, `s62_sample_size_chr22.json`,
`s63_plink_export_chr22.json`, `s64_saturation_chr22.json`, `s65_identity_ablation_chr22.json`,
`s66_control_provenance_chr22.json`,
`refine-logs/EXPERIMENT_PLAN_DEEP_20261016.md`,
`refine-logs/EXPERIMENT_RESULTS_DEEP_20261016.md`

## Overall Verdict: PASS

## Integrity Status: pass

The auditor replayed all ten S61 cells in memory (all 120 membership AUROCs
reproduced within 1e-7; all ten GWAS rows within 1e-8), re-parsed both cohorts from
the local VCF (genotypes, labels, positions matched), verified disjoint roles for
all ten splits, and replayed both prior-round synthetic witnesses. No project files
were modified. It found no fabricated ground truth, no output-derived normalization
that inflates a score, and no honest-private-matrix access in the adaptive ranking.

Three rounds of findings were raised and closed:

- Round 1 (verdict FAIL): the GWAS row fitted on a combined honest+non-member panel
  (fixed to the honest study alone); a `sqrt(recovered_rank)` divisor on the
  membership score (removed); a missing conditioned-AF baseline and permutation null
  (both added); target directions selected by shadow-PC index (replaced by
  residual-against-passive-span); heterozygosity measured on a separate swap
  (measured on the same releases); the results summary's marker/honest-n/passive-flat
  statements (corrected).
- Round 2 (verdict FAIL): the conditioned baseline was not composition-balanced
  (replaced by equal-quota removal per group); default-`pinv` reconstruction
  instability (single relative cutoff now drives the basis, the pseudo-inverse and
  every derived metric); adaptive release accounting (output 1 is now a shared
  passive release, so q costs q releases for every strategy); conditioned-AF gaps
  absent from the JSON (added).
- Round 3 (verdict WARN): adaptive-minus-random CI endpoints printed in the summary
  (added; all include zero); lambda upper-bound rounding (1.28); the adaptive
  ranking's default-`pinv` prior span (replaced by the same relative cutoff).
- Final round (verdict PASS): the focused real-data replay reproduced all twelve
  checkpoint AUROCs and ranks exactly; continuous-metric differences at most
  2.73e-8; every printed number reconciles with the JSON.
- Round-2 re-audit (`01a14455-c02b-7c93-b103-24ed8aace24a`, verdict FAIL → WARN →
  **PASS**): the expanded artifacts (chr1, s62, s63, s64, s65, s66) were replayed
  fresh. Closed findings: a chr1 3pop q4 "5/5" contradiction (now 3/5); the
  adaptive-minus-random "excludes zero only at ..." list (now three cells); a
  mis-attributed adaptive maximum (now chr22 EUR q2); an unsupported "clipping
  changes only the random-span control" (now all three arms); a missing s64 tracker
  row (added); the reviewer-only raw-eigenvector ablation (now artifact s65); an
  S66 coordinate-system bug (span-fitting and scoring now share one transform); an
  overstated "requires the identity" (now scoped to the tested raw-V attack); an
  overstated "all at chance" (now the descriptive 0.455–0.509 bound); a mislabeled
  s65 baseline field; and two rounding slips. Final verdict PASS **within the
  declared experimental scope** — it does not establish broader generalization or
  preprocessing equivalence.

## Checks

### A. Ground Truth Provenance: PASS
Membership labels are real 1000G study-inclusion labels on a disjoint 5-role split;
the private truth $C_H$ is computed by the evaluator separately from the attack
measurements. The panel phenotype is simulated and used only by the GWAS-power row.
`code/s61_spectral_tomography.py:242,253,384`.

### B. Score Normalization: PASS
The membership score is a plain norm `np.linalg.norm(X @ B, axis=1)`; no divisor
from the recovered rank remains (`code/s61_spectral_tomography.py:275`). The other
denominators (visible_share by `||C_H||_F`, target errors by `|x^T C_H x|`) are
reference-based.

### C. Result File Existence: PASS
The seven result JSONs exist (s61 chr22/chr1, s62, s63, s64, s65, s66); every number in
the results summary reconciles with them (322 table values checked in the final round-2
pass, all matching).

### D. Dead Code Detection: PASS (one unused library helper)
`metrics`, `witness`, `gwaspower`, `af_baselines`, `summarize`, the permutation
null, the bootstrap helper and the adaptive ranking are all called and their outputs
appear in the JSON. `common.stratified_permutation_auroc` is unused (s61 implements its
own permutation test) — unused library functionality, not a fabricated metric.

### E. Scope Assessment: WARN (absent claim language; open comparisons named)
The report uses no "comprehensive/robust/extensive". The un-run planned comparisons
(per-release scientific-cost curve, full federated pipeline) are named as not
supported. Two chromosome prefixes × five seeds; one k=10, block length 20, candidate
budget 24; q=16 only in the chr22 saturation cell; chr replication reuses subjects and
is not independent-cohort replication.

### F. Evaluation Type: real_gt (membership), simulation_only (phenotype row)
Membership AUCs and coverage metrics use real-genotype experimental inclusion
labels; the GWAS-power row is a simulated phenotype on real genotype covariates;
the witness checks are self-supervised consistency checks.

## Central mathematical claim
$Y_t=V_t\Lambda_t-K_tV_t=C_HV_t$ is algebraically correct and `Y` is constructed
from the released eigenpairs, not from `C_H @ V`. The three metrics follow distinct
paths; the adaptive ranking reads only partner data, the shadow panel, past released
eigenvectors and the public center/scale.

## Claim Impact
- D1 (identity $Y=C_H V$ holds on real layout, incl. PLINK export read-back): supported.
- D2 (AF-preserving swap keeps marginal statistics exact): supported.
- D3 (coverage/recovery grow with #outputs; passive flat): supported.
- D4 (accumulated measurement beats the strongest passive baseline): supported.
- D5 (honest-analysis summary statistics stable across releases): supported (summary-statistic scope).
- D6 (adaptive beats random): not supported (max +0.014, ~+0.010 or less elsewhere).
- D7 (scatter recovery survives export precision): not supported (float64-only).
- D8 (the increment is specific to the tomography identity): supported (s65).
- D9 (the no-release control gives no separation across preprocessing choices): supported (s66).
