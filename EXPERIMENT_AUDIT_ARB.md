# Experiment Audit Report — Round 22 (fair comparison & interface confirmation)

**Date**: 2026-10-17
**Auditor**: Codex / gpt-6-astra, read-only (cross-family; thread `01a14948-947a-7f22-9862-5dd2d9935d3c`)
**Project**: /data/lvvyongji1/GWAS
**Artifacts audited**: `code/s67_joint_passive.py`, `code/s68_interface.py`,
`code/s69_reliable_dirs.py`, `code/s70_known_fraction.py`, `code/s71_protection_cost.py`,
and their `results/raw/s67_*_chr{22,1}.json`, `s68_interface_chr22.json`,
`s69_*_chr{22,1}.json`, `s70_*_chr{22,1}.json`, `s71_*_chr{22,1}.json`;
`refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md`,
`CLAIMS_FROM_RESULTS_ARB_20261017.md`, `refine-logs/EXPERIMENT_PLAN_ARB_20261017.md`.

## Overall Verdict: PASS (for the narrowed claim)

## Integrity Status: pass

The auditor independently reproduced every persisted s71 summary (within 9.6e-15),
all 30 random-q8 oracle comparisons, all 24 s70 fraction summaries, and the chr22 EUR
s69 precision results; it also replayed s71 with only the tested markers and reconciled
the corrected lambda-GC / false-positive-rate values.

## Checks

### A. Ground Truth Provenance: PASS
Membership labels come from assigned inclusion roles; the scatter reference is the
honest genotype. The s71 effect reference was corrected from a distribution SD to the
sampled per-marker effect, and the effect measure is now a scale-free Pearson
agreement.

### B. Score Normalization: PASS (after fix)
s69 previously derived a reported rank from `sv / sv[0]` on the joint operand; it now
measures the reliability floor from the export perturbation of the recovered operand
`Y_t = V_t Λ_t − K_t V_t` and no reported AUROC/TPR is self-normalized.

### C. Result Existence and number match: PASS
All JSON pairs exist; the s69 stale-schema file (joint operand, unreliable-rank-0) was
replaced. The auditor's independent reconstruction matched every summary exactly after
the reporting corrections.

### D. Dead Code: WARN
`common.stratified_permutation_auroc` remains unused. s70's over-budget `fused_cal`
diagnostic was removed from the execution path (it trained on all calibration labels,
exceeding the known-member budget).

### E. Scope: WARN
Two chromosomes × five seeds; one k=10; chr1 reuses the same 1000G subjects (not an
independent cohort). The interface manifest demonstrates LOCAL software capability under
an assumed access model, now labeled as such.

### F. Evaluation Type
real_gt (membership) / synthetic_proxy (scientific phenotype) / self_supervised_proxy
(scatter geometry).

## Action Items
None remaining; all round-22 findings (s67 baseline selection, s69 operand, s70
fraction denominator and balanced removal, s71 marker universe / effect reference /
recovered-count table, scope wording, and the center-claim narrowing) are closed.

## Claim Impact
- J1 (float64 oracle-baseline superiority): supported.
- J2 (not center/scale; not fixed-input resampling): supported.
- J3 (interface items demonstrated; access assumed): supported as capability.
- J4 (approx solver agreement): supported as a lower bound.
- J5 (precision-reliable directions): supported as a heuristic bound.
- J6 (increment across known-member fraction): supported.
- J7 (scientific cost): partial.
- J8 (independent subjects): not run.
- J9 (bounded-precision conjunction in printed output): not established.

## Verdict history
FAIL → FAIL → FAIL → FAIL → WARN → **PASS** (six review packets on one thread).
