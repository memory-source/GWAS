# Round-16 Review — limitations pass (run_20261008_25d5b534)

**Reviewer**: pi (artifact-level verification). Independently recomputed the size-resolved
crossing brackets from `s21_sensitivity.json`; verified the three paper edits, the control
ranges, the CLAIMS_FROM_RESULTS update, and the compile.
**Date**: 2026-10-08

## The three changed items — verified

1. **Size-resolved crossing brackets** ✓. Paper: "[0.027,0.052] at n=1600, [0.052,0.077] at
   n=800, [0.103,0.151] at n=300 (s21_sensitivity.py)". My independent recomputation
   (last-positive / first-negative measured fst_like per size) gives exactly
   [0.0267, 0.0518] / [0.0524, 0.0770] / [0.1029, 0.1506] — an exact match, and "the crossing
   value of structure rises as n shrinks" is correct. The anchoring claim is sound: the
   scale-matched real cells (PC wins at Fst 0.0997, n=60 — below the n=300 bracket's lower
   bound) and the unconfounded 3/5-pop cells (gaps ≈0/−0.03 at n=450/715, near the interpolated
   crossing) are consistent with the simulator's surface, and the text says exactly that
   ("the real cohorts sample that surface at small n") without overclaiming.
2. **chr2 3-population cell** ✓. Now "a mild PC>AF exception: |gap| ≤ 0.023 at the crossing,
   where the chromosome-level sign is unstable---chr1 is −0.006 there" — chr1 3pop is −0.0055 ✓.
   The false "small size" explanation is gone (paper and results doc both now say "proximity to
   the crossing"). This is the correct explanation, and it matches the bracket framing (the
   3-pop rung sits near the crossing, where the sign is operating-point sensitive).
3. **Eigenvalue-only range "0.55–0.62"** ✓ — exact for the within-population controls across
   s49 (0.545–0.596) and s50 chr1/chr2 (0.563–0.618, covering chr1:YRI 0.618). The
   reference-basis range "0.55–0.60" (s49 0.552–0.595) and oos "0.55–0.58" (0.551–0.581) are
   also now exact. All round-15 range imprecisions resolved.

Additional regression checks: "n_study=28–36" corrected in paper and results doc (B1/B4) ✓;
intro now reads "monotone in measured structure **under a size-controlled F_ST sweep**" with the
stale "below the real-data level" clause removed ✓; "both sides"/"measured-structure law": 0
hits ✓; CLAIMS_FROM_RESULTS.md updated to the round-13 verdict (C6 yes-scoped / C7 yes / C8 yes)
with an honest "what the results do not support" section ✓; compile 17 pp, 0 errors, 0
undefined refs ✓. No remaining text-vs-JSON inconsistency found in this pass.

## 1. Score: 8.5 / 10

## 2. Limitations a top-venue reviewer would raise (ranked; blocker vs disclosed scope bound)

**L1 (largest; disclosed scope bound — not a blocker).** No real biobank-scale real-genome
cohort exists in the paper: the size-resolved law's large-n/low-structure quadrant is carried
by the parametric simulator, the real large-n side is the mixed 3/5-pop panels (n=450/715)
plus the semi-synthetic mosaic cohort, and real within-population cohorts are capped at
≈30–36 members by 1000G population sizes. Whether AF beats PC at Fst≈0 *given* large n is
undetermined on real data (CLAIMS_FROM_RESULTS discloses exactly this). Inherent to available
data; the paper labels the semi-synthetic arm and scopes the real claims to the sampled
surface. Not a blocker.

**L2 (moderate; partially disclosed — the one cheap strengthening left).** The real round-13
gaps carry no per-split uncertainty (5-seed means only; the JSONs aggregate). The near-crossing
cells (3pop −0.0055/+0.0231, chr22 −0.005/−0.030) have sign instability the text now
acknowledges qualitatively ("at the crossing … sign unstable"), but a referee can ask for
per-split values or a bootstrap CI. Fix if desired: emit per-split rows in the s48/s50 JSONs
(one-line change to `_agg`) — not a blocker, since every affected claim is hedged.

**L3 (minor; disclosed).** The PC-arm orientation on the real ladder uses evaluation labels
(`_cal` = max(A,1−A)), inflating near-chance values by ~±0.016; headline cells (≈1.0, 0.65–0.69)
are unaffected. Disclosed in CLAIMS_FROM_RESULTS. Not a blocker.

**L4 (minor; disclosed).** The real-data M axis is one thinning contrast (s30) plus the
per-arm densities; the s48/s49/s50 arms run at one density each. Not a blocker.

**L5 (minor; disclosed).** f-DP ε_lb is resolution-limited — ordering-only wherever the
empirical ROC hits a zero-FN corner (valid flags in the JSON; rule-3 floor; explicit
"ordering information" wording). Not a blocker.

**L6 (minor; disclosed).** Governance is a plaintext release-object prototype with a stated TCB
(local user with process control; parent runtime); the native-code boundary is now tested
(3/3 ctypes denial), not asserted. Not a blocker.

**L7 (minor; disclosed).** The missingness channel is constructed-by-design (ρ = member-heavy
batch fraction; ρ=0 null exactly at chance 0.502). Not a blocker.

**L8 (nit).** 3 seeds/cell in s37/s41/s42 (tight stds on the headline deltas). Not a blocker.

## 3. Text-vs-artifact inconsistencies

None found. The three changed items match their artifacts exactly; the previously flagged
imprecisions (n_study 28–36; control ranges; the chr2 explanation) are all corrected; no stale
"law"/"both sides"/monotone-for-real-ladder wording remains; CLAIMS_FROM_RESULTS matches the
round-13 artifacts.

## 4. Minimum fixes

None required. Optional: L2's per-split dispersion (one-line `_agg` change + re-run or a
bootstrap over the existing caches) if a referee asks.

## Verdict

**READY (yes).** The paper now states exactly what its artifacts show: the joint size×structure
gating anchored to the size-resolved simulator brackets, real genomes on both the within-pop
(PC>AF, mechanism verified, controls at chance) and across-pop (AF-first at large n) sides, the
chr2 exception correctly attributed to the crossing, and every quoted number matching its JSON
key. The remaining limitations are disclosed scope bounds inherent to the available data, not
defects. 8.5 is the honest ceiling given L1; nothing found this pass blocks submission.
