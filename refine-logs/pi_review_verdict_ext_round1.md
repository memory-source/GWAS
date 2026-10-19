# Review verdict — extension round 1 (s51/s52/s53 chromosome extension)

**Reviewer**: pi (independent, artifacts read directly; framing not trusted)
**Scope**: s51 chr2, s52 chr1/chr2, s53 chr1 additions + paper/claims/tracker consistency
**Verdict: READY — 8.5/10**

The extension does what it claims, every quoted number pairs 1:1 to a live JSON
key, the previously broken planted control now recovers, and no headline claim
changes. Remaining weaknesses are minor (internal-doc wording, one completeness
gap, one optional-coverage item).

## Direct answers to the four attention items

1. **chr2 wording vs numbers — no overclaim in the paper.** The paper says chr2
   "reaches near-parity at the top rung (+0.109/+0.059/+0.018/+0.006)" — all four
   rungs are positive in `s51_real_lowstruct_chr2.json` (0.1087/0.0593/0.0176/0.0057),
   and "near-parity" is the correct description. No "crossing" is claimed for chr2
   anywhere in main.tex. The crossing bracket "n_study≈260–400" is carried by
   chr1 (~258, between rungs 2–3) and chr22 (~385, between rungs 3–4); chr2 does
   not cross within the ladder and the paper does not say it does. One internal
   doc (tracker L1) does say "ladder crosses: … chr2 …" — wrong for chr2, see W1.
2. **s52 planted control recovers on all three chromosomes.** Real arm at chance
   (chr22 .488/.496, chr1 .497/.497, chr2 .497/.498); planted 1500-marker copy
   detected (chr22 .798/.797, chr1 .776/.788, chr2 .772/.779) — sitting at the
   60%-members-frac construction ceiling (~0.80), i.e. every planted member is
   detected and unplanted members stay at chance. The round-15 fix is in the code
   (si-only planting via `plant_phased_ibd(H0[si],…)`; donor is a study member and
   stays in the released set; non-members scored against the same release, no
   contamination; `seg=1500`, `thr=0.78` shared across chromosomes — no
   per-chromosome retuning). Control certified, so the real-arm negative is
   meaningful.
3. **Paper-number ↔ JSON audit: all match.** Verified 1:1: s51 ladder values and
   Fst (chr22 +0.120/+0.042/+0.007/−0.026 at Fst 0.0911→0.0873 = "0.087–0.091";
   chr1 +0.058/+0.002/−0.010/−0.018 at Fst≈0.10; chr2 as above at Fst≈0.086);
   s51 pooled gaps (chr22 +0.075/+0.051/+0.044; chr1 +0.110/+0.165/+0.048;
   chr2 +0.181/+0.111/+0.154) and the n_study=167–220 cap; the crossings bracket;
   s52 all 12 quoted AUROCs + thr 0.78 + 1500-marker; s53 "0.48–0.52" (fixed
   bases 0.484–0.517, incl. random 0.478–0.523) and full loadings 0.808–0.885.
   Spot-checked prior-round ranges the abstract leans on: s48 twelve single-pop
   cohorts +0.051..+0.136 (n_study 28–36), 2-pop +0.114..+0.198, 3pop +0.007@450 /
   5pop −0.041@715, eigonly 0.376–0.597 ("0.38–0.60"), s50 within-pop
   +0.062..+0.091. No mismatch found.
4. **Extension changes no headline claim.** C1/C2/C6/C7 scoping is unchanged and
   strengthened: the size ladder is now 3-chromosome (chr22/chr1 cross, chr2
   near-parity), the haplotype negative is 3-chromosome with a working control,
   and the fixed-eigenvalue-basis negative covers chr1. The new largest real
   PC>AF cells (s51 pooled chr2, up to +0.181) sit in the small-n in-sample
   regime (n_study 167–220) already disclosed by the mechanism section; the
   abstract's "+0.05 to +0.14" range is explicitly scoped to the twelve s48
   cohorts and remains correct. Cross-script agreement is good (s53 chr22/chr1
   EUR full-loadings 0.829/0.867 vs s51 pooled 0.820/0.864 — two independent
   code paths, consistent).

## Integrity checks run (beyond number pairing)

- s52: all 60 per-split AUROCs lie exactly on the 240×240 midrank lattice
  (2·AUC·240² ∈ ℤ); aggregates equal per-split means to 1e-12; embedded params
  (per=80, thr=0.78, seg=1500) identical across the three chromosome JSONs and
  match the current script constants.
- s51: gap = pc_auc − af_auc exact in all 21 cells; pooled n_pool 503/504/661 =
  true 1000G phase-3 super-population counts; ladder n_study 120/240/360/480;
  held-out-half orientation via shared cal_mask (calibration rows excluded from
  evaluation) confirmed in `s48.cohort_cell`/`_orient_auroc`.
- Plant mechanics verified in `s32_plant_phased_ibd`: donor = study member
  (idx[0] per super-pop), planted rows only in si, leave-self-out removal indices
  correct for the vstack layout, donor hap stays in the release (the channel the
  control certifies).
- s53 lam_weight control uses the study's own released eigenvalues (1/λ_j on the
  first k standard axes) — the divide-by-own-statistic is the control under test
  and lands at chance, with full loadings positive in the same cells (0.81–0.89):
  the harness detects signal when present.
- mtime anomalies benign: s52 script edited 23:50 between chr22 (23:47) and chr1
  (23:51) runs and s53 script at 02:58 after its 02:34 JSON — embedded JSON
  parameters match current constants in both cases (docstring-level edits).
- Docs propagated: EXPERIMENT_RESULTS L1/L3/L4, TRACKER L1/L3/L4, and CLAIMS C6/
  C16/C17 all carry the chr2/chr1 extension numbers; root EXPERIMENT_AUDIT.md
  refreshed (02:59) and its three text nits are verifiably applied in main.tex.

## Remaining weaknesses (ranked; none blocking)

- **W1 (low, must-fix-in-docs): tracker L1 wording "ladder crosses: … chr2
  +0.109/+0.059/+0.018/+0.006".** chr2 never crosses (top rung +0.006). The
  paper and CLAIMS are clean; the tracker is not. MINIMUM FIX: reword the L1
  note to "ladder falls monotonically (chr22/chr1 cross; chr2 near-parity at the
  top rung)". Optionally append "on chr22/chr1" after "crossing locations
  (n≈260–400)" in main.tex §scale, EXPERIMENT_RESULTS L1, and CLAIMS C6 (chr2's
  crossing, if any, is an extrapolation beyond 480, not an interpolation).
- **W2 (low): seeds paragraph completeness.** The Operating-point parenthetical
  lists split counts for the older arms but not s51/s52/s53; "real-data …
  mechanism studies 20 splits" could be misread as covering s53 (5 splits).
  MINIMUM FIX: add "the fixed-basis probe, the second haplotype statistics, and
  the real size ladder 5 splits/seeds per cell" to the parenthetical. (s51's 5
  splits is stated in §scale; s52/s53 counts appear nowhere in the paper.)
- **W3 (low, residual statistical honesty): the ladder story rests on 5-split
  means with per-rung sd 0.02–0.09; the near-zero rungs (chr22 per120 +0.007±
  0.026, chr2 per160 +0.006±0.023) are individually indistinguishable from zero,
  so the crossing locations carry wide uncertainty.** The paper already hedges
  ("interpolations of the five-split means", "near-parity") and the direction is
  consistent across three chromosomes and with the simulated boundary, so this
  is a caveat, not an error. MINIMUM FIX (optional): quote the per-rung sd for
  the ladder gaps once (values already in the JSONs) or state that the top rungs
  are within split noise of zero — one sentence.
- **W4 (optional hardening): s53 real cells cover only EUR/EAS on chr22/chr1**
  (no AFR, no chr2). The paper's scoping ("on simulated and real chr22/chr1") is
  accurate, so this is not an error. If a reviewer pushes, the extension is a
  one-line `chroms=(22,1,2)` + AFR loop and a cheap re-run.
- **W5 (note, no action): near-identical operating points straddle zero across
  arms** — s9 chr2 mixed cohort −0.063 at (n=450, M=4664) vs s51 chr2 per160
  +0.006 at (n=480, M=3084). Both differences (fewer markers, slightly larger n)
  push the gap in the direction the paper's density/size framework predicts, and
  s50 already labels chr2's 3-pop cell as sign-unstable at the crossing; no
  contradiction, but a one-line footnote reconciling M across arms would pre-empt
  a careful reviewer.

## What plainly holds

- The chr2 extension is real, consistent, and correctly worded ("near-parity",
  not "crossing").
- The s52 planted control now recovers everywhere at its construction ceiling
  with a shared threshold — the round-15 blocker is genuinely closed.
- The s53 eigenvalue-weighted fixed basis is a well-designed control (positive
  reference in the same cells) and its "0.48–0.52" quote is exact.
- No paper number mismatches its JSON; internal-doc propagation is complete
  except W1.
- The extension is strengthening-only: it adds external validity to C6/C7/C16/
  C17 without disturbing any scoped claim.

**Score: 8.5/10 — READY.** (Prior closure round: 8/10. The extension adds
verified multi-chromosome coverage and a working positive control; the half-point
is withheld for W1's internal-doc wording and W3's unquantified rung uncertainty,
both one-sentence fixes.)
