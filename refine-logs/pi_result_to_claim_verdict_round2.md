# Result-to-Claim Evaluation — Round 2 (after closing the round-1 evidence gaps)

Date: 2026-10-04
Evaluator: pi (paper-writing skill, claim-strength + experiment-gate rules)
Method: every cited number re-read from `results/raw/*.json` (not trusted from framing);
calibration/split protocols inspected in `code/s9_orientation.py`, `s10_paired.py`,
`s11_robust.py`, `s12_crossover.py`, `s13_real_loadings.py`; carried-over numbers
spot-checked against `s1c_frontier.json`, `s2c_mechanism.json`, `s1d_bias.json`,
`s5_boundary.json`, `s7_real.json`, `s8_families.json`.

---

## 1. claim_supported: **partial**

The claim is now a *scoped* partial — much stronger than round 1. Five of six
components are individually well-supported with the new evidence, and every
round-1 gap was addressed. What keeps the verdict at partial is that the new
experiments themselves revealed a genuine empirical boundary: **the headline
comparative ("PC diagnostic is THE dominant channel") is contradicted in ordering
on the only real dataset tested, and is regime-dependent even in simulation**
(S12 crossover). As stated — "In multi-site GWAS, the released study-fitted PC
diagnostic is the dominant intermediate membership-leakage channel" — the claim
is over-broad. The scoped version (below) is fully supported.

## 2. what_results_support

Per sub-claim:

**(a) PC > AF among aggregate releases, ~0.88 vs ~0.84 — SUPPORTED in simulation.**
- S10 (`s10_paired.json`): paired Δ(PC−AF) = +0.0461, 95% CI [0.0372, 0.0550],
  t p≈4.5e-10, Wilcoxon/sign p≈5.96e-8 (= 2⁻²⁴, the discrete minimum for n=25 —
  all 25 seeds agree in sign). Calibration verified clean in code: orientation
  chosen on a train half only; combined logistic fit on train, scored on test;
  no test-side label use. This is the formal paired test round 1 lacked, and it
  is sound.
- Convergent: S2c pcdist_only 0.8813 vs F5_full 0.8446 (paired 5/5 seeds);
  S8 pc_distance 0.873–0.889 vs f1_af_lrt 0.836–0.859 (paired 5/5 per seed).

**(b) "+0.10 increment" — SUPPORTED in simulation at the studied operating point.**
- S10: combined vs AF = +0.1030, 95% CI [0.0965, 0.1099], p≈0.
- But see (3) below: operating-point-specific.

**(c) Mechanism = in-sample-ness of the study-fitted spectrum — SUPPORTED (strongest new evidence).**
- S9d (`s11_robust.json`): a site's own fitted basis separates its members at
  0.8813/0.8800; the *other* site's basis scores 0.5099/0.5047 ≈ swapped-null
  (0.5047/0.5099). Out-of-sample transfer collapses to chance. This is a
  decisive transfer test, not just a correlation.
- Converging: S2c T2 independent-spectrum control Δ+0.003≈0 vs T3
  released-spectrum control Δ+0.10 preserved; permuted labels 0.49;
  pcalign-only 0.50.
- Caveat: all mechanism evidence is simulation; the real-data analog was not run.

**(d) Individual leakage is real, not an ancestry confound — SUPPORTED.**
- Disjoint-negative controls at chance everywhere (0.5019–0.5047 sim; 0.5036 real).
- S9b proves the control is not a no-op: the reference-basis control sits at
  0.5019 on matched members/non-members but rises to 0.6341 when the basis
  separates true ancestry groups — so chance-level control is informative, not broken.
- Real data: S7 per-population attack 0.808 vs structure control 0.504.

**(e) Released-object granularity governs the leak — SUPPORTED on BOTH sim and real (best-reproduced component).**
- Sim (S9c): full loadings 0.8816; top-k monotone (0.527/0.607/0.655/0.740 for
  k=1/2/3/5); eigenvalue-only 0.5100 ≈ chance.
- Real (S13, 20 splits): full loadings calibrated 0.6522; monotone in k
  (0.509/0.526/0.545/0.579); eigenvalue-only 0.5145 ≈ chance; AF-only 0.6879.
  The ordering reproduces qualitatively on real data — this closes the
  "granularity" gap at both levels.

**(f) Self-reference bias ≈ +0.13 — SUPPORTED but narrow.**
- S1d: 0.1319 mean (per-seed 0.119–0.143, 5 seeds, simulation, F5_full tier only).
- Unchanged from round 1; still sim-only and single-tier.

**Gap-closure check (round-1 list → status):**
1. Real-data PC orientation study + matched-operating-point sim → S9 ✓
2. Formal paired significance test PC vs AF → S10 ✓ (sim-only)
3. Crossover map → S12 ✓
4. Explicit release-tier taxonomy → closed definitionally in the paper
   (`paper/main.tex` defines "intermediate" and the F-tier ladder) + the
   empirical granularity ladder in S9c/S13
5. Batch-correlated missingness → S9a ✓ (0.7614 ± 0.0116 — a real, moderate
   channel: above chance, below PC and AF at this operating point; fills the M3 gap)

## 3. what_results_dont_support

1. **The unscoped dominance claim.** On real 1000G chr22 (S9, 10 disjoint
   3-way splits), sign-calibrated PC = 0.6533 < AF-only = 0.6921 — orientation
   "closer" on 10/10 splits (verified `orientation_closer_frac=1.0`), so this is
   not an orientation artifact. S13 confirms: full-loadings calibrated 0.6522 <
   AF-only 0.6879 on 20/20 splits. The PC channel does not dominate on the only
   real dataset tested. "Dominant in multi-site GWAS" as a general statement is
   falsified by the paper's own new artifacts.
2. **The magnitudes as general numbers.** ~0.88 vs ~0.84 are simulation numbers.
   Real-data magnitudes are 0.65 vs 0.69 — lower and ordered the other way.
3. **"+0.10" as a general increment.** S12 shows the PC−AF increment spans
   +0.149 (400×500) to −0.006 (1600×8000) across the (panel × marker) grid — a
   13-fold-plus range that crosses zero. The +0.10 (and +0.103 CI) holds at the
   standard operating point only. On real data the increment of AF+PC over AF is
   +0.043 (0.7348 vs 0.6921) — positive in 10/10 splits but untested formally.
4. **S12 cell precision.** 6 seeds per cell, means only, no stds/CIs in the
   grid — the crossover pattern is systematic but per-cell uncertainty is
   unreported.
5. **Framing corrections (do not change conclusions, but fix before writing up):**
   - S9 real arm has **10** split seeds, not 20 (`meta.real_split_seeds=10`;
     script default `tuple(range(10))`). "Orientation closer on 20/20" is
     10/10 for real; 20/20 applies to the sim arm.
   - S12 peak is +0.1493 at **400×500** (`n400_m500`), not at 400×1000
     (+0.1477 there). Trivially different but the stated location is wrong.
   - S9's "oracle" is attack-side `max(naive, flipped)` **using labels**
     (documented in the script); on real data `pc_oracle == pc_flipped` exactly
     and `naive + flipped = 1.0` (sign-flip symmetry). Report it as attack-side
     sign calibration, not as an upper bound requiring no labels.
   - Cross-script sign-convention inconsistency: S9-sim `pc_naive` is
     anti-predictive (0.09) while S10's `train_pc_naive` is predictive (0.88) —
     the two scripts define the naive orientation oppositely. Internally
     consistent, but any reader cross-checking S9 vs S10 will see contradictory
     "naive" AUROCs; document the convention.
   - S10's config (1200/group, 4000 SNPs) is the same regime as but not
     identical to S2c/S8 (1000/group, 4000 SNPs); af_only 0.8324 vs F5_full
     0.8446. Fine, but cite the config.

## 4. missing_evidence

1. **Real-data paired CI for the AF+PC increment** (+0.043, 10/10 splits positive,
   no formal test). Cheap: a paired t/bootstrap on S9's per-seed values.
2. **Mechanism on real data.** S9d is sim-only; an S9d analog on chr22 (split the
   study into disjoint sites, in- vs out-of-sample basis) would close the
   mechanism-on-real gap.
3. **Second real dataset.** All real-data conclusions rest on 1000G chr22 with
   1363 SNPs. The PC<AF ordering could be chr22- or thinning-specific. chr1 or
   HGDP replication needed before any real-data generalization.
4. **Non-parametric simulation for the PC channel.** S4 (msprime) only tested the
   AF tier (0.642). PC dominance under a coalescent generative model is unknown;
   all PC evidence is Balding-Nichols + injected LD.
5. **S12 uncertainty.** Stds/CIs per grid cell, more seeds at the crossover cells.
6. **M3-batch sweep.** Batch-correlated missingness evaluated at one operating
   point (0.7614); its ranking vs AF/PC across the S12 grid is unknown.

## 5. suggested_claim_revision

The paper already scopes dominance to the simulated panels (`main.tex`:35–36,
84–86) — keep that scoping and add the new boundaries. Revised claim:

> "In simulation (Balding-Nichols + injected LD, 3-way disjoint panels), the
> released study-fitted PC diagnostic is the largest single aggregate leakage
> channel at our standard operating point: 0.88 AUROC vs 0.83 for the AF summary
> release (paired Δ+0.046, 95% CI [0.037, 0.055], n=25, p<10⁻⁸); combining
> channels yields +0.103 over AF (CI [0.097, 0.110]). This dominance is
> regime-dependent: the PC advantage peaks in the small-panel/mid-marker regime
> (+0.149 at 400 individuals × 500–1000 markers) and vanishes when panels are
> large and marker density high (Δ−0.006 at 1600×8000, where AF saturates).
> On real 1000G chr22 data the PC diagnostic alone does not dominate
> (sign-calibrated 0.653 < AF 0.692, 10/10 splits), but adds a small consistent
> increment on top of AF (0.735 vs 0.692), and the granularity ordering — full
> loadings > top-k loadings > eigenvalues ≈ chance — reproduces on real data
> (0.652 > 0.509–0.579 ≈ 0.51). Mechanism: in-sample-ness of the study-fitted
> spectrum (own basis 0.88, other site's basis 0.51 ≈ null; released-spectrum
> control preserves the increment, independent-spectrum control kills it).
> Self-reference scoring inflates naive estimates by ≈0.13 in simulation."

Key moves: (i) "dominant" → "largest single aggregate channel at the studied
operating point, regime-dependent"; (ii) state the real-data ordering flip
honestly; (iii) anchor +0.10 with its CI and operating point; (iv) promote the
granularity result to the best-reproduced component (it holds on sim AND real);
(v) keep "+0.13" scoped to simulation.

## 6. next_experiments_needed

Priority order:
1. **Real-data held-out-site mechanism test** (S9d analog on chr22) — cheapest,
   closes the biggest remaining mechanism gap.
2. **Second real dataset** (chr1 or HGDP) for the PC<AF ordering — determines
   whether the flip is chr22-specific.
3. **Paired significance for real-data AF+PC vs AF** (+0.043) — trivial compute.
4. **msprime cross-check for the PC channel** (S4 infrastructure exists; run the
   PC attack under the coalescent model).
5. **S12 with stds/CIs** and more seeds at 1600×4000/8000 cells.
6. **M3-batch sweep** across the S12 grid.

## 7. confidence: **high**

All 14 new cited values verified against raw files; all calibration and split
protocols inspected in code (S10's held-out-half calibration confirmed clean;
S9's oracle confirmed label-using and documented); carried-over numbers
spot-checked against round-1 artifacts. The verdict rests on convergent,
pre-declared, controlled evidence with formal paired statistics — the residual
uncertainty is about generalization (one real dataset, parametric sim), not
about whether the numbers are real.

---

## Evidence table (all verified)

| Claim component | Key artifact | Value | Status |
|---|---|---|---|
| PC > AF (sim, paired) | s10_paired.json | Δ+0.0461, CI [0.0372, 0.0550], 25/25 sign | SUPPORTED (sim) |
| +0.10 increment (sim) | s10_paired.json | +0.1030, CI [0.0965, 0.1099] | SUPPORTED (sim, op-point) |
| Regime dependence | s12_crossover.json | +0.149 (400×500) → −0.006 (1600×8000) | SUPPORTED (sim) |
| PC < AF (real) | s9_orientation.json, s13_real_loadings.json | 0.6533 < 0.6921; 0.6522 < 0.6879 | CONTRADICTS unscoped claim |
| Mechanism (in-sample-ness) | s11_robust.json s9d | 0.8813 own vs 0.5099 other ≈ null | SUPPORTED (sim) |
| Not ancestry confound | s11_robust.json s9b, s7_real.json | control 0.5019 (matched) vs 0.6341 (structure-present); real 0.808 vs 0.504 | SUPPORTED |
| Granularity governs leak | s9c (sim) + s13 (real) | sim 0.882/0.51; real 0.652/0.515, monotone in k | SUPPORTED (both) |
| Self-reference +0.13 | s1d_bias.json | 0.1319 (5 seeds, F5_full) | SUPPORTED (sim, narrow) |
| Batch-correlated M3 | s11_robust.json s9a | 0.7614 ± 0.0116 | New channel filled (one op-point) |
| AF-tier real ordering | s7_real.json | perpop 0.808 > full 0.684 > 1dp 0.564 | unchanged, consistent |
