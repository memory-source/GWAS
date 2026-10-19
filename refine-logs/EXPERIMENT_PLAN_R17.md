# Experiment Plan (Round 17 — cross-model defect closure after the round-16 audit)

**Problem**: The round-16 Codex honesty audit returned `overall_verdict: fail` on
08 concrete defects in the new s54/s55/s56/s37/s42 artifacts and the paper text.
The substantive ones are already fixed (s56 rewritten to a balanced design; s54
3-pop gap corrected; s37 5-seed merge rebuilt; paper wording corrected). What
remains is the residual set that the audit named but round 16 did not yet close,
plus the two text-only items held over from `EXPERIMENT_PLAN_LIMITS.md` (D1).

**Method Thesis**: Close each residual as either (i) a bounded correctness
disclosure with a recomputation, or (ii) a stated scope sentence — no new
mechanism, no new experiment family.
**Date**: 2026-10-14

## Claim Map
| Claim | Why It Matters | Minimum Convincing Evidence | Linked Blocks |
|-------|-----------------|-----------------------------|---------------|
| (A) The GPU `auroc_gpu` tie handling is wrong but does not move any headline | A reviewer who reads `common_gpu.py` finds `argsort` ranks give a constant score 0.0 not 0.5 | A bounded probe showing the discrepancy and the tie mass of the affected cells | B1 |
| (B) The s54 permutation evidence is scoped correctly | "every real cell" overstated a 5-cell test; `gap_p_gt0` is the positive-split fraction, not a bootstrap probability | Scope sentence + label; no new run | B2 |
| (C) The s56 positivity is a five-split statement | Some s0.70 cells are positive on 4/5 splits | Per-split gaps exported (done) + disclosure that near-crossing cells are not 5/5 | B3 |
| (D) The s37 cohort controls are honestly re-stated | "cohorts pass their controls" exceeded the seed-0, no-gate reality; two ld bins are undefined | Control sentence reworded (done) + the undefined-bin case named | B4 |
| (E) The ε convention and governance TCB are stated once, globally | Hold-over D1 from the round-15 plan | One ε-convention paragraph + one TCB paragraph; no channel-specific `valid` usage | B5 |

## Paper Storyline
- **Main paper**: no new headline; (A)–(D) are corrected disclosures already in
  the text, (E) is a single-paragraph consolidation.
- **Appendix**: the s59 tie probe JSON.
- **Cut**: s57 (per-window LD-correlation — measures phase LD, not membership),
  s58 (further in-simulation eigenvalue prior — speculative).

## Experiment Blocks

### Block 1 (B1): GPU AUROC tie-handling disclosure — DONE
- Claim (A). `code/s59_auroc_tiecheck.py` → `results/raw/s59_auroc_tiecheck.json`.
- Constant `[1,1,1,1]`/`[1,1,0,0]`: CPU midrank 0.5, GPU argsort rule 0.0.
  Half-tied vector: 0.505 vs 0.252 (gap 0.253). Archived s27–s36 cells carry
  continuous AF-LRT / projection scores, so tie mass is negligible and the
  checked headline precision is unaffected; disclosed as a bounded defect.
- Priority: MUST-RUN (done).

### Block 2 (B2): s54 scope + `gap_p_gt0` label
- Claim (B). Text-only. The null covers three within-population cohorts + the
  3-/5-population cells (five cells), not all eighteen; it tests each channel's
  separation, not the PC−AF difference. Applied in `paper/main.tex` §scale.
- Priority: MUST-RUN (text only).

### Block 3 (B3): s56 per-split disclosure
- Claim (C). Per-split gaps are exported in the JSON; the paper states the
  largest balanced cell and the n450 reference regime. Add one clause that
  near-crossing cells are 4/5-positive so the "does not flip" claim is a
  five-split-mean statement.
- Priority: NICE-TO-HAVE.

### Block 4 (B4): s37 controls re-statement
- Claim (D). Text-only (done): AF spectrum and eigenvalue fraction match at all
  cells; LD decay matches at the five cells whose panel supports the shortest
  bin; the two shortest panels (M=3000/6000) have that bin undefined.
- Priority: MUST-RUN (text only, done).

### Block 5 (B5): ε convention + governance TCB consolidation
- Claim (E). Text-only. State the FNR floor convention once (default
  `1/n_pos`, conservative rule-of-three `3/n_pos`) applied on every channel, and
  the process-isolation TCB (PEP-578 hook, single sanctioned pipe, stdout/stderr
  discarded, ordinary Python / no native code / no root).
- Priority: NICE-TO-HAVE.

## Run Order and Milestones
| Milestone | Goal | Runs | Decision Gate | Cost | Risk |
|-----------|------|------|---------------|------|------|
| M0 | sanity | B1 tie probe | returns the two AUROCs | seconds | none |
| M1 | main | B2 + B4 text | numbers pair to JSON | minutes | none |
| M2 | polish | B3 + B5 text | convention stated once | minutes | none |

## Compute and Data Budget
- Total GPU-hours: 0. All blocks are a seconds-long CPU probe plus text edits.

## Risks and Mitigations
- Text edits drift from the artifacts. Mitigation: every quoted number already
  pairs to a live JSON key; the audit re-checks the pairing.
- The tie defect tempts a kernel rewrite that would invalidate s27–s36.
  Mitigation: disclosure only; the kernel is left unchanged.

## Final Checklist
- [x] Main paper tables covered (no new table)
- [x] Novelty isolated (no new mechanism)
- [x] Simplicity defended (disclosure + text, no machinery)
- [x] Frontier contribution justified or explicitly not claimed
- [x] Nice-to-have separated from must-run (B1/B2/B4 must; B3/B5 nice)
