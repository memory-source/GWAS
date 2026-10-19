# Code review verdict — round 2 (verify fixes)

Scope: verify the two round-1 fixes (CRITICAL no-op ancestry control; MAJOR maf dropping
f>0.475) in `/data/lvvyongji1/GWAS/code/` against the fresh JSONs in
`/data/lvvyongji1/GWAS/results/raw/` (all re-run 2026-10-02 17:39–17:54, matching the code
fix timestamps 17:38–17:53).

## Verdict: PASS

Both fixes are correct, both old defects are gone from the active pipeline, the structure-only
control is meaningful, and the paper's 0.502-at-chance claim matches the fresh data. No remaining
CRITICAL/MAJOR issue would invalidate the numbers. Four MINOR text/hygiene items below.

---

## Check 1 — both fixes correct, both old defects gone: CONFIRMED

**maf() fix.** `common.py` now `f = G.mean(0)/2.0; return np.minimum(f, 1 - f)` with a docstring
documenting the old bug; `qc_mask` therefore keeps f ∈ [0.05, 0.95]. Both inline sites in
`s4_msprime.py` fixed identically (downsample + attack mask). Effect verified in fresh JSONs:
S0 `n_snps_tested` = **3814/4000** (was 1874; the request text says 3817 — 3814 is what the JSON
shows, immaterial), `n_causal_tested` = **9/10** (was 4). S4 now tests 2326–2466 of its 4000
downsampled variants (was ~half).

**Ancestry-control fix.** New `anc_only_attack(target_G, study_loadings, study_mean, mask)`:
centers the target by the passed (public) study mean, projects onto the passed loadings, scores
−‖projection‖. Panel/label alignment verified correct in both callers:
- `s0_sanity.py`: now runs the **three-way disjoint split** (300/300/300 per ancestry,
  n_per_group=1000) and attacks **nonmem** (third panel) as negatives; the control scores
  members (label-1 block) and nonmem (label-0 block) with the study-fitted basis/mean —
  correct panels, correct labels, correct ordering.
- `s1c_frontier.py`: two controls per tier, both with correct panel/label alignment —
  `ancestry_only_auroc` = **reference-fitted basis** (load_ref/mean_ref; site and nonmem both
  out-of-sample) and `study_pc_distance_auroc` = **study-fitted basis** (load_s/mean_s; the S2c
  leak, naive orientation).

**Old defects gone.**
- The no-op helper `ancestry_only_attack` is marked deprecated and is **not called by any active
  script** (s0/s1c/s2c/s4 import `anc_only_attack`). It survives only in the stale comparison
  artifacts `s1_frontier.py` / `s1b_disjoint.py`, whose JSONs are not cited as current results
  (see MINOR-2 for the one abstract claim that leans on them).
- The dead `(n,1,m)` line in `af_lrt_attack` is removed (only the correct broadcast remains).
- `s2c_mechanism.py` docstring now says "FARTHER" — verified correct: members' in-sample
  projections sit farther from the centroid (mean pcdist more negative than non-members'),
  matching the paper's "independent panel projects with systematically smaller spread".

## Check 2 — structure-only control meaningful; paper's 0.502 claim correct: CONFIRMED

Fresh `s1c_frontier.json`: `ancestry_only_auroc` = **0.5020 at every tier** (matches Table 1
caption "at chance (0.502) at every tier"); per-ancestor F5 = 0.826/0.810/0.799 (caption ✓).

The control is meaningful, not vacuously chance: both target panels are projected onto the
**independent** reference panel's PC basis and centered by the reference mean, so both are
out-of-sample — no study-fitting leak can enter by construction. Verified it has genuine (if
modest) power when structure exists: scoring an ancestry-mismatched pair (mixed-ancestry site vs
single-stratum non-members) gives 0.487/0.483 vs matched 0.500–0.521 across seeds — the control
responds to real structure, it is not a constant.

Frontier and mechanism numbers all match the fresh JSONs: S1c F2 0.584 / F4 0.604 / F3 0.788 /
F5 0.795 / F5-per-site 0.935 (TPR@1% 0.024/0.032/0.108/0.121/0.398, CIs match); S2c base 0.839,
pcdist-only 0.881, pcalign 0.501, base+M2 0.941 (Δ+0.102), T2 Δ+0.003, T3 Δ+0.102, eig 0.9998,
perm 0.486; S2b base 0.844, +PC 0.942 (Δ+0.098); S4 0.642 full / 0.649 round-1; S3b parity ✓;
S5/S6 re-run. The abstract's numbers (0.88 / 0.79–0.80 / 0.94 / 0.94; ≈0.2 self-reference
inflation) all trace to these files. The ≈0.2 self-reference claim rests on the old-vs-disjoint
comparison (s1_frontier 0.916 vs s1b_disjoint 0.720 on the same simulator config), which was
same-maf in both arms, so the Δ is self-consistent.

## Check 3 — remaining issues: no CRITICAL/MAJOR; four MINOR

1. **MINOR (paper text, §controlled "Attack" paragraph).** It still says "plus an ancestry-only
   control (distance to the **study** PC centroid) as a confound upper bound" — stale. The
   control actually reported in Table 1 is the **reference-basis** control (the caption describes
   it correctly). Fix: update the Attack paragraph to "projected onto an independent reference
   panel's PC basis", and adjust "confound upper bound" (it is a structure-only null; the
   study-fitted distance is the leak, already reported in §mechanism).
2. **MINOR (stale artifacts).** `s1_frontier.py` / `s1b_disjoint.py` still call the deprecated
   no-op and ran with the old maf. The abstract's "≈0.2" self-reference claim is from those
   old-maf runs; the comparison is internally consistent (same maf both arms) but re-running it
   with the fixed maf would make it airtight.
3. **MINOR (S0 single-seed lambda).** Fresh S0 `lambda_gc_uncorrected` = 1.17 (corrected 0.95) —
   the plan's S0 claim holds directionally in this realization, but the number is
   realization-unstable: across seeds λ_unc ranges **1.2–6.6** in both old (n_per_group=700) and
   new (1000) configs (verified). The old 6.35 was realization noise, **not** the maf bug (the
   buggy maf gives 1.07 on the same fresh population). Consider reporting the median over seeds
   or noting the variance.
4. **MINOR (S0 JSON naming).** s0's `ancestry_only_auroc` key (= 0.115) is the **study-fitted**
   PC distance at naive orientation — i.e. the S2c leak (oriented ≈0.885), not an
   ancestry-confound control. Rename to `study_pc_distance_auroc` for consistency with s1c. Not
   cited in the paper.

Also noted (no action): the two s1c controls are tier-invariant (they use the full QC mask, not
the tier's released AF), so "at every tier" repeats one number — true, but the caption could say
"tier-invariant control".

---

## Bottom line

- **PASS.** maf() fixed correctly everywhere (3814/4000 tested; 9/10 causal); `anc_only_attack`
  correctly projects onto study-fitted loadings with correct panels/labels in S0 and S1c; the
  deprecated no-op is out of the active pipeline; the dead line and docstring are fixed.
- The structure-only control is meaningful (out-of-sample by construction, responsive to real
  structure), and the paper's 0.502-at-every-tier claim matches the fresh JSON exactly.
- All fresh JSON numbers match the paper's abstract, Table 1, Table 2 and §mechanism.
- Apply the four MINOR text/hygiene fixes at convenience; none invalidates the reported numbers.
