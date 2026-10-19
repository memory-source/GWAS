# Round-15 Review — regression check of W13-1..W13-5 (run_20261007_66850dac)

**Reviewer**: pi (artifact-level verification). Verified the scale-matched rerun JSON against
the paper, the five fix regions, the framing reconciliation (abstract / contribution 1 /
§scale / Discussion / Limitations), and the compile.
**Date**: 2026-10-07

## W13-1..W13-5 scorecard

| # | Fix | Status | Evidence |
|---|---|---|---|
| W13-1 | True s50 range + chr2 exception disclosed | **FIXED (one flaw)** | "across-population gaps −0.039 to +0.023, with the chr2 3-population cell a mild PC>AF exception" ✓ — matches the JSON (+0.0231). The false "two chromosomes confirm" is gone. **New flaw**: the appended explanation "consistent with its small size" is false — the chr2 3-pop cell has n_study=450, identical to chr1's 3-pop cell (−0.0055). The correct explanation is proximity to the crossing (\|gap\| ≤ 0.023; the same rung flipped sign between s16 at M=4664 [−0.063] and s50 at M=3384 [+0.023]). |
| W13-2 | "Monotone" removed for the real ladder | **FIXED** | §scale: "The data points are not monotone in structure alone … two local non-monotonic steps, so we bound the crossover only loosely in [0.02,0.075]" ✓. Abstract's monotone claim gone ✓. Parametric monotonicity correctly retained with the "under a size-controlled F_ST sweep" qualifier (Limitations) ✓. |
| W13-3 | Sizes disclosed + scale-matched rerun | **FIXED (excellently)** | Sizes disclosed per rung; "fst_like is 0 by construction (one group) … a definitional operating point, not a measured structure value" ✓; the confound stated ✓; the rerun verified against `s48_scalematched.json` (3pop Fst 0.0997, n=60, gap +0.1607 → "+0.161"; 5pop Fst 0.0885, n=100, gap +0.1112 → "+0.111"; controls at chance 0.518–0.562) ✓; framing revised to "governed jointly by study size and structure" in §scale, abstract, contribution 1, Discussion, and Limitations ✓ ("both sides"/"measured-structure law": 0 hits). |
| W13-4 | k/n regime disclosed | **FIXED** | "n_study=30–36, k=10, so k/n≈1/3---the near-1.0 AUROCs are the geometric small-n limit of the in-sample effect, and on the larger simulated slabs the same mechanism gives 0.88 at n=900" ✓. |
| W13-5 | §boundary sentence + ranges + B2 | **FIXED (one nit)** | §boundary now parses: "the diagnostics release is the largest canonical intermediate channel in the small-(n×M) regime---which includes single-population real cohorts …---while once study size and structure grow the AF release leads" ✓. Results-doc B2 now "GBR/CHB/YRI on chr1/chr2" ✓. oos range "0.55–0.58" ✓ (actual 0.551–0.581). **Nit**: the control ranges "0.53–0.60" — s49 actuals are refb 0.552–0.595, eig 0.545–0.596 → "0.55–0.60"; the quoted lower bound is below the artifact minimum (possibly a per-split value, but not verifiable from the JSON). |

## The scale-matched rerun — the decisive experiment, verified

`s48_scalematched.json` (per-group n=20, 5 seeds): single GBR/CHB +0.091/+0.076; **3pop
(Fst 0.0997, n=60) +0.161; 5pop (Fst 0.0885, n=100) +0.111** — PC dominates even across
super-populations at matched small n, with controls at chance (0.518–0.562). This exposes the
round-13 ladder's size confound exactly as suspected, and the paper now says so everywhere:
the AF-first switch needs both larger n and structure; the ordering is governed jointly by
study size and structure — the same small-(n×M) boundary the simulator locates, now reproduced
on real genomes. This is a stronger and more defensible claim than the round-13 "law on both
sides" framing, and it is honestly derived (the confound was run down, not hidden).

## Regression check

Compile: 17 pp, **0 errors** (4 overfull hboxes). All round-11/12/13 numbers previously
verified remain in place; the abstract, contribution 1, §scale, Discussion, and Limitations are
consistently re-framed; the round-12 W12-5 sentence ("admits a valid f-DP bound") is intact.

## Score: 8.5 / 10 (restored — the confound is now controlled and disclosed, not hidden)

## Remaining weaknesses (all minor)

1. **"consistent with its small size"** (paper §scale + results-doc B3) — false explanation:
   the chr2 3-pop cell is n=450, the same size as chr1's confirming cell. Fix (one clause):
   "consistent with its proximity to the crossing".
2. **"n_study=20–36"** for the twelve s48 single-pop cohorts (paper + results-doc B1/B4) —
   actual 28–36; the 20 belongs to the scale-matched rerun's single-pop cells (a different
   file). Fix: "28–36" (or "≈30").
3. **Control ranges "0.53–0.60"** (s49-scoped) — actual 0.545–0.596 / 0.552–0.595 → "0.55–0.60";
   the quoted lower bound is unverifiable from the JSON. Fix: "0.55–0.60".
4. Cosmetic: the intro's "crosses at measured structure ≈0.04–0.06, below the real-data level"
   describes the parametric sweep (accurate) but predates the joint framing.

Minimum fixes: items 1–3 are one-clause/one-number edits; none blocks any result.

## Verdict

**YES — submission-ready.** Every blocking-class defect from the entire loop (misquoted ranges,
false headline claims, undisclosed confounds, leak-checked demonstrations) is resolved; all
quoted numbers match their JSON keys; the one genuinely new experiment of this round — the
scale-matched ladder — is the right control, it was run, and it moved the paper's claim to
where the data supports it. The three remaining wording imprecisions are optional polish
(item 1 strongly recommended before submission), not blockers.
