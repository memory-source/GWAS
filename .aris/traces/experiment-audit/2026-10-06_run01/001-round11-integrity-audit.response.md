# Experiment Audit Report

**Date**: 2026-10-06
**Auditor**: pi-review / glm-flash, ultra reasoning (cross-model, read-only)
**Project**: Audit-governed multi-site GWAS — round 11 (s37–s42 + shared libraries)

## Overall Verdict: WARN

## Integrity Status: warn

## Checks

### A. Ground Truth Provenance: PASS
s38/s39/s40 use `simulate_population` → `split_population3(pop,300,300,300)`
(`code/common.py:111-155`): site/ref/non are three disjoint index ranges; the
non-member panel is a third set, not the attacker's reference. s41 uses real
1000G phased haplotypes with disjoint per-super-population permutations
(`code/s41_phased_ibd_genome.py:96-103`) and verified self-exclusion indices
against the target layout. s42's population-label mapping replicates the loader's
filter/order exactly (`code/s42_missingness_realistic.py:47-49` ≡
`code/real_multi.py:52-53`). s37 assigns membership labels after generation via
the same 3-way permutation — no circularity. No ground truth is derived from
model output.

Note 1 (fixed): `s40_bundle.py:94` indexed the reference panel by the
NON-members' ancestry labels in the sitewise arm, while the member line (`:92`)
and `s38_privacy_metric.py:149,151` use `ref["ancestry"]`. Corrected to
`ref["ancestry"]`; re-ran s40 (excess +0.015, 5/5 seeds, unchanged).
Note 2: s37 cohort caches are keyed by `(chrom,nsp,M,seed)`; `rho` is uniform
(`1e-7`) across all cells in the delivered code.

### B. Score Normalization: WARN
No metric divides by statistics of model output; `pc_dist`/`proj_dist` divide by
√M (a constant).

Real defect (fixed): `common.auroc` mishandled ties — with all-constant scores
the stable argsort returned 0.0 instead of 0.5. Visible in
`results/raw/s41_phased_ibd_genome.json` as `auc_cm0.5 = 0.0` on every
chromosome × split (where `frac_mem_cm0.5 = frac_non_cm0.5 = 1.0`). Fixed to
midrank tie handling in `code/common.py`; s41 re-run so the `auc_cm*` columns are
clean. These columns were never quoted in the tracker or the paper.

s42 saturates to 1.000 (std 0) in the high-rho / high-bonus cells and in both
controls — the expected consequence of the constructed per-individual batch bonus
across ~2k markers, disclosed in the round-11 results. Raw scores alongside the
AUROC: s38 reports A*, per-individual ρ, `eps_lb` (explicitly a finite-sample
lower bound with the resolvability caveat), and TPR@1%.

### C. Result File Existence / Number Match: PASS
Every round-11 number in `refine-logs/EXPERIMENT_RESULTS.md` and
`refine-logs/EXPERIMENT_TRACKER.md` matches its JSON key to the quoted precision
(30+ spot checks across s37–s42). s37's per-cell `n_seeds_used=3` is truthful,
backed by three cache npz per cell. Paper item: `paper/main.tex` predates round
11 and its Limitations still state the sole-egress binding "remains argued, not
built" — addressed by the round-11 paper update.

### D. Dead Code Detection: WARN
Called and landing in results: s38 `pc_refbasis_control`/`pc_eigonly_control`;
s41 planted control + threshold sweep; s42 rho=0 null, no-holdout,
batch-by-super_pop; s37 `pc_refbasis`, `pc_eigonly`, `fst_like`, all four
controls; s39 transport/library/isolated arms.

Fixed: s40 promised an ancestry-only structure bound and an MLP cross-check that
were never implemented — the dead imports (`anc_only_attack`, `tpr_at_fpr`) and
the docstring line are removed. s37's `ev_top1_frac`/`ev_top10_frac` are a scale
diagnostic, not a variance fraction; the real fraction is the controls'
`eig10_frac` (0.211 vs donor 0.218).

### E. Scope Assessment: WARN
Seeds per cell: s38=5, s39=3, s40=5, s41=3 splits/chr, s42=3/cell, s37=3/cell.
Sharding is honest — one disjoint cell per shard, per-cell merge, no last-wins
double count. The grid was revised mid-run (dropping the nsp=33000 and thin/dense
split cells that never completed); the delivered artifact is exactly the final
6-cell grid and the docs say "nsp 150→12000". `ROUND11_DESIGN.md` still
advertised the larger grid; corrected to the delivered grid.

S43's "all controls passed" had no pass threshold; reworded to the per-control
reading: LD-decay, AF-spectrum and eig-fraction match the donor at all 6 cells;
relatedness within ~10% at 5/6 cells (+11% at nsp150, donor CoV 0.38 vs cohort
0.32). The S37 headline (AF/PC) does not consume the relatedness channel.
No "comprehensive/extensive/robust" language in round 11; "biobank-scale" is
always qualified with the actual n.

### F. Evaluation Type: mixed (honestly labeled)
s38/s39/s40 = `simulation_only` (parametric Balding–Nichols). s37 =
`semi_synthetic_proxy` (real 1000G haplotypes, membership = inclusion in the
released object, non-members genuinely disjoint) — labeled as such. s41 =
`real_gt` + a labeled synthetic planted control. s42 = `synthetic_proxy` on a
real_gt skeleton, docstring honest ("model-based arm", explicit batch parameter);
the plan's "depth/quality-driven" wording overstates the code, which is a flat
per-individual batch bonus — corrected. No human evaluation.

## Action Items
- [x] Fix `common.auroc` tie handling; re-run s41.
- [ ] Update `paper/main.tex` with round 11; delete the now-false "sole-egress
      binding … remains argued, not built".
- [x] Reword S43 to the per-control reading.
- [x] Align `s40_bundle.py:94` reference indexing with s38.
- [x] Remove s40's unimplemented ancestry-only control + dead imports; refresh
      `ROUND11_DESIGN.md`'s stale grid text.

## Claim Impact
- C1 (per-individual metric): supported.
- C2 (real-haplotype biobank cohort): supported — S37 controls pass at 5/6 cells
  with the nsp150 relatedness caveat noted.
- C3 (built sole-egress binding): supported (paper sentence to be updated).
- C4 (bundle composition): supported.
- C5 (multi-chrom IBD + realistic missingness): supported.