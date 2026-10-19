# Experiment Plan (Round 12 — residual-limitation closure)

**Problem**: Round 11 closed five limitations, but the round-11 result-to-claim
judge (`CLAIMS_FROM_RESULTS.md`, verdict `yes`) and the integrity audit list four
remaining, explicitly-scoped gaps: (a) the f-DP `eps_lb` on the four high-AUROC
channels is floor-driven (`valid=False`, no confidence-corrected FNR floor); (b)
no per-individual metric (A*/rho/eps) on real or semi-synthetic data — every C1
number is simulated; (c) the largest delivered biobank cell is nsp=12000, not the
nsp=33000 biobank corner the design targeted; (d) the bundle excess is one
simulated manifest with an oracle-selected single object, and the C3 sole-egress
binding is a scoped demonstration with no native-code boundary test.
**Method Thesis**: Each residual is a *measurement-bound* statement, not a new
mechanism — so each is closed by re-running an existing estimator at a new
operating point or against a second configuration, not by new machinery.
**Date**: 2026-10-07

## Claim Map
| Claim | Why It Matters | Minimum Convincing Evidence | Linked Blocks |
|-------|-----------------|-----------------------------|---------------|
| C1b The per-individual metric is not a simulation artefact | C1 currently holds only on simulated panels; a privacy reviewer wants it on real structure | A*/rho/eps_lb on real 1000G (chr22) and on the s37 semi-synthetic cohort, ordering preserved, controls at ~0 | B1 |
| C1c The f-DP `eps_lb` is quotable, not floor-driven | `valid=False` on the top channels undermines the absolute number | A confidence-corrected FNR floor (rule-of-three) so `valid=True` at some delta on the top channel | B1 |
| C2b The AF-first ordering survives at the biobank corner | "nsp=12000" leaves the extrapolation open | One cell at the largest feasible nsp (≤24 GB GPU) still AF-first | B2 |
| C4b The bundle excess is not one manifest / an oracle pick | One manifest + oracle max is a single-point result | A second manifest (different objects/precisions) and an uncoarsened per-site AF single channel | B3 |
| C3b The sole-egress binding's stated boundary is tested | The threat model names "no native code"; untested, it is asserted | A native-code escape attempt (ctypes dlopen) is blocked / the boundary is confirmed | B4 |

## Paper Storyline
- **Main paper must prove**: C1b and C1c (they strengthen the C1 headline directly).
- **Appendix can support**: C2b, C4b, C3b (each tightens an existing scoped caveat).
- **Experiments intentionally cut**: a real biobank / UKB panel (not obtainable);
  an MPC/SPU implementation (out of scope — only the egress boundary is built); a
  full 22-autosome IBD sweep (chr1/2/22 already bound it; diminishing returns).

## Experiment Blocks

### Block 1 (B1): Per-individual metric on real + semi-synthetic data + confidence floor
- Claim tested: C1b, C1c.
- Why it exists: round-11 C1 is simulation-only and its eps_lb is floor-driven.
- Dataset / split / task: (i) real 1000G chr22 (existing `real_multi.get_real_panel_chr`,
  the same 3-way split as `s16`); (ii) an s37 semi-synthetic cohort at nsp=2000
  (real haplotype structure). Re-express the AF-LRT and study-fitted-PC channels.
- Compared systems: the same estimators as `s38_privacy_metric.py`
  (`membership_advantage`, `per_individual_posterior`, `eps_lower_bound`), with a
  new `fnr_floor="rule3"` option: `FNR_floor = 3/n_pos` (the rule-of-three
  upper confidence bound on a zero observed FN count) instead of `1/n_pos`.
- Metrics: `(AUC, TPR@1%, A*, eps_lb(delta) with rule-3 floor, validity flag,
  median rho, LR_max)`, per channel, plus the two controls.
- Setup: CPU (`env/bin/python`); reuse `prior`/`alt` classes of `s38` verbatim.
- Success criterion: ordering preserved on real and semi-synthetic data (PC<AF on
  real as in `s16`; AF leads on the cohort); controls at A*≈0; `valid=True` at
  delta=0.05 on the top channel with the rule-3 floor.
- Failure interpretation: if rule-3 still gives `valid=False` at every delta, the
  eps_lb is reported as ordering-only everywhere (an honest, tighter bound).
- Table / figure target: extend the round-11 per-individual table; one appendix row each.
- Priority: MUST-RUN.

### Block 2 (B2): Largest feasible biobank cell
- Claim tested: C2b.
- Why it exists: the delivered grid stops at nsp=12000.
- Dataset / split / task: one s37 cell at nsp=16000 (measured feasible: 44 s,
  fits the 24 GB GPU with the subset-fitted PC), M=3000, 3 seeds.
- Compared systems: the `s27`/`s37` battery (af_only, pc_calibrated, delta,
  pc_refbasis, pc_eigonly, in-sample gap, fst_like) + the four controls.
- Metrics: the s27 cell keys.
- Setup: GPU (`/opt/miniforge/bin/python3`); `gram_topk_subset` at n_fit=20000.
- Success criterion: AF still leads (delta<0) and controls pass at nsp=16000.
- Failure interpretation: if PC overtakes AF at nsp=16000, the boundary is
  re-stated to exclude the biobank corner (a real revision).
- Table / figure target: one added row in the round-11 cohort table.
- Priority: MUST-RUN.

### Block 3 (B3): Second bundle manifest + uncoarsened single channel
- Claim tested: C4b.
- Why it exists: one manifest + an oracle-selected single object is a single point.
- Dataset / split / task: `release_objects` under a second, plausible manifest
  (per-site AF at full precision + eigenvalues + lambda_GC@2dp + the call-rate
  table — omit pooled AF and p-values); and a variant whose strongest single
  object is the *uncoarsened* per-site AF (the oracle max at full precision).
- Compared systems: the joint out-of-fold logistic bundle vs max single, with the
  same structural-zero control and leave-one-out attribution.
- Metrics: `bundle_agg_auc`, `max_single_auc`+argmax, `excess`, `loo_drop`, noise control.
- Setup: CPU, minutes.
- Success criterion: the excess is positive in 5/5 seeds under the second manifest
  and under uncoarsened-single too.
- Failure interpretation: if the excess vanishes when the single channel is
  uncoarsened, the round-11 excess is a coarsening artefact — reported as such.
- Table / figure target: extend the round-11 bundle table.
- Priority: MUST-RUN.

### Block 4 (B4): Sole-egress native-code boundary test
- Claim tested: C3b.
- Why it exists: the threat model names "no native code"; the PEP-578 hook is
  advisory against native code and this is currently asserted.
- Dataset / split / task: extend `enforce.py`'s injection set with a `ctypes`
  escape attempt (`ctypes.CDLL(None).write(...)` / `dlopen`) issued from the child
  study process.
- Compared systems: the isolated runtime's audit-hook outcome on the native-code
  route vs the library prototype.
- Metrics: `native_egress_blocked` (bool), plus the existing bypass matrix.
- Setup: CPU, seconds.
- Success criterion: the audit hook raises on `ctypes.dlopen`/`dlsym` (the route is
  denied) — OR, if a route escapes, that escape is reported and the threat model is
  narrowed to "no ctypes/native loading" explicitly.
- Failure interpretation: an escape is a real finding; the paper's threat-model
  sentence is tightened to match.
- Table / figure target: one row in the round-11 bypass matrix.
- Priority: MUST-RUN.

## Run Order and Milestones
| Milestone | Goal | Runs | Decision Gate | Cost | Risk |
|-----------|------|------|---------------|------|------|
| M0 | metric + enforcement sanity | B1 rule-3 unit check, B4 native test | rule-3 lowers eps at a fixed ROC; ctypes route resolves | <5 min CPU | native route API mismatch |
| M1 | per-individual on real + cohort | B1 full | ordering preserved; controls at 0 | ~15 min CPU | real split differs from s16 |
| M2 | bundle second manifest | B3 full | excess>0 both manifests | ~10 min CPU | second manifest under-covers |
| M3 | biobank corner cell | B2 (nsp16000 ×3) | AF leads; controls pass | ~5 min GPU | GPU contention |
| M4 | write-up + re-audit | all | every number traces to a JSON | — | none |

## Compute and Data Budget
- Total GPU: <0.5 h (only B2). Everything else CPU (minutes).
- Data: existing 1000G caches (chr22) + one new s37 cohort cache (nsp=2000, 16000).
- Human evaluation: none.
- Biggest bottleneck: the nsp=16000 cohort generation + PC fit (subset-fitted).

## Risks and Mitigations
- rule-3 floor still `valid=False`: report eps_lb as ordering-only everywhere; the claim narrows honestly.
- Second manifest under-covers: if the bundle collapses, report the (non-composition) negative — a legitimate result.
- Native-code route unavailable on this host: fall back to documenting the boundary and marking the test "not run" (do not claim a block that did not happen).
- GPU contention on nsp=16000: shard lightly; the cell is 44 s when idle.

## Final Checklist
- [ ] Per-individual metric on real + semi-synthetic data (C1b)
- [ ] Confidence-corrected eps_lb (C1c)
- [ ] Biobank-corner cell (C2b)
- [ ] Second manifest + uncoarsened single (C4b)
- [ ] Native-code boundary test (C3b)
- [ ] Nice-to-have separated from must-run
