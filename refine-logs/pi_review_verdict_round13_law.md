# Round-13 Review (run_20261007_66850dac)

**Reviewer**: pi (artifact-level verification). Read s48/s49/s50 code + all three JSONs, the
four paper edit regions, the results/plan write-ups, and the reconciliation state of the older
real-data sections.
**Date**: 2026-10-07

## The four adversarial questions — answered

**(2) Within-population PC-dominance: genuine mechanism, but the "law" framing overreaches.**
The controls are real and at chance: `cohort_cell` fits the reference-basis control on the
disjoint 30-person ref panel (both scored panels out-of-sample for it) and the eigenvalue-only
control on a random orthonormal basis; across all 16 s48 cohorts refb=0.514–0.587 and
eigonly=0.511–0.601, chance-consistent at each rung's own size (at n=450/715: 0.508–0.524,
within ~1 SE). So the within-population signal is NOT within-population structure — it is the
in-sample-ness of the study-fitted basis, exactly the paper's claimed mechanism. BUT: (a) the
within-population PC AUROCs are 0.987–1.000 at n_study≈30 with k=10 (k/n≈1/3); s49's
samehalf/oos arms fit ~15 samples with k=10 (k/n≈2/3) — the ≈1.0 is the geometric extreme of
the in-sample mechanism at tiny n (the paper's own prior own-basis results: 0.88@900 sim,
0.685@450 real). The regime is never disclosed. (b) `fst_like` for single-pop cohorts is 0 by
construction (one group ⇒ hs=ht) — the §scale text says "between-population Fst≈0", which is
accurate, but it is definitional, not a measurement of within-population structure.

**(3) s49 held-out-half: no train/test leak — the demonstration is real.** oos arm: basis fit on
`site[half]` (~15), scored targets `site[~half]` (the OTHER half, not in the fit) vs non-members
→ 0.551–0.581. samehalf arm: basis fit on the scored half itself (intended in-sample condition)
→ 0.998–1.000. Same scored members, only the basis differs — the design cleanly isolates
in-sample-ness. Genuine, but the k/n≈2/3 regime must be stated for the ≈1.0 to be credible.

**(1) Number check:** s48 (12 single-pop gaps +0.050…+0.134 ✓; the four quoted AF/PC pairs
match exactly ✓; 2-pop +0.112…+0.172 ✓; 3/5-pop −0.005/−0.030 at Fst 0.075–0.087 ✓); s49
(top-3 0.881–0.959 ✓, full 0.994–1.000 ✓, eig 0.545–0.596, refb 0.552–0.595, oos 0.551–0.581,
samehalf 0.998–1.000 ✓); s50 within-pop +0.051…+0.086 all positive ✓. Compile 16 pp, 0 errors ✓.

**(4) Overclaims found — three, one serious:**
- **The quoted s50 range hides a contradicting data point.** §scale: "across-population gaps
  −0.005 to −0.039". The JSON: chr1 3pop −0.0055, chr1 5pop −0.0380, chr2 5pop −0.0391, and
  **chr2 3pop +0.0231 (PC > AF at Fst 0.0855, n=450)**. The +0.023 point is outside the quoted
  range and contradicts "across super-populations AF leads"; "Two independent chromosomes
  confirm it" is false as stated (chr2's 3-pop rung does not confirm).
- **"The gap is monotone in measured Fst" (§scale) / "with the gap monotone in measured
  structure" (abstract) — contradicted by the artifact**: single-pop (Fst=0) gaps +0.05..+0.13 →
  2-pop (Fst 0.004–0.020) gaps +0.11..+0.17 (gap RISES with Fst); 5pop (Fst 0.0747) −0.030 →
  3pop (Fst 0.0874) −0.005 (gap RISES with Fst). Two non-monotonicities. The plan's own success
  criterion ("gap monotone-decreasing in measured Fst") was not met.
- **"Crossing near Fst≈0.05" is an interpolation** — no rung exists between Fst 0.020 and 0.075;
  the honest statement is that the sign flips between 0.02 and 0.075.
- **The structure-vs-size confound is unacknowledged**: the ladder's panel sizes co-vary with
  structure (n_study 28–36 → 52–66 → 450 → 715; JSON `n_study`). The within-pop rungs cannot be
  enlarged (1000G population sizes), but the decisive size-matched rung — a 3-super-population
  cohort at per=20 (n≈60, Fst≈0.09), a one-parameter change to `cohort_cell` — was not run.
  Without it, the real ladder shows "both orderings occur on real genomes" (true, valuable) but
  cannot isolate structure as the governing variable; the isolation currently rests on the
  size-controlled parametric sweep (s21), which the text does not say.
- **§boundary "scope the claim" sentence is mangled**: "We therefore scope the claim: the
  diagnostics release is the largest canonical intermediate channel as a measured-structure
  law---…" — the old scoping sentence and the new law were spliced into a claim that no longer
  parses (the diagnostics release is not "a measured-structure law").
- **Control/oos range imprecisions**: §scale "controls sit at 0.53–0.60" (s49 actual 0.545–0.596
  → "0.55–0.60"); "transfers at chance (0.55–0.59)" (actual 0.551–0.581 → "0.55–0.58");
  "eigenvalue-only 0.51–0.60 in every real cohort" — s50 chr1 YRI eigonly = 0.618 exceeds it
  (chance-consistent at n=36, but the quoted range is wrong). Results doc B2 says "GBR/CHB
  chr22" for s49 — the delivered s49 ran chr1/chr2 only (JSON chroms=[1,2]).

## What genuinely holds (and is valuable)

Real genomes now exhibit both orderings, in the directions the parametric law predicts: 12/12
single-population and 5/5 two-population cohorts show PC>AF; 3 of 4 across-population cells show
AF≥PC; the controls (independent basis, random basis) are at chance, so the within-population
signal is the released basis' in-sample-ness — the paper's own mechanism — and s49 shows it
needs the loadings (top-3 0.88–0.96 → full ≈1.0, eig-only at chance) and does not transfer
(0.55) while the same-half basis gives 1.0. This is a real strengthening of the paper's story.

## Score: 8.0 / 10 (from 8.5: the new evidence is valuable, the law framing outruns the artifact)

## Remaining weaknesses (ranked) and minimum fixes

**W13-1 (major — misquoted range / hidden exception).** "across-population gaps −0.005 to
−0.039" omits s50 chr2 3pop = +0.0231. **Fix**: quote the true range (−0.039 to +0.023), state
the chr2 3-pop exception explicitly, and replace "Two independent chromosomes confirm it" with
"chr1 replicates both sides; chr2 replicates the within-population side, with its 3-population
rung near the crossing (+0.023)".

**W13-2 (major — false monotonicity + interpolated crossing).** Abstract "monotone in measured
structure"; §scale "monotone in measured Fst, crossing near Fst≈0.05". **Fix**: replace with
"positive at measured Fst ≤ 0.02, negative at Fst 0.075–0.09, with the crossover bracketed in
[0.02, 0.075]"; drop "monotone" for the real ladder (keep it for the size-controlled parametric
sweep, where it holds); keep the plan's failure-interpretation honesty.

**W13-3 (major — confound unacknowledged).** The ladder co-varies size (28–36/52–66/450/715)
with structure. **Fix**: disclose the sizes in §scale and state that within-population cohorts
are capped at ≈30 members; run the one-parameter size-matched control (3-pop cohort, per=20,
n≈60 — `cohort_cell(G, lab, per=20)`), report it either way, and scope the real-data claim as
"consistent with the measured-structure law, whose size-controlled isolation is the parametric
sweep" unless the matched rung pins the crossing.

**W13-4 (moderate — regime disclosure).** State n_study≈30, k=10 (k/n≈1/3) for s48 and the
~15-sample fits (k/n≈2/3) behind s49's samehalf/oos arms, and that the ≈1.0 is the in-sample
mechanism at its small-n extreme (prior own-basis magnitudes: 0.88@900, 0.685@450).

**W13-5 (minor).** Fix the mangled §boundary sentence ("the diagnostics release is the largest
canonical intermediate channel as a measured-structure law" — no longer parses); correct the
control ranges ("0.55–0.60" for s49's controls; "0.55–0.58" for oos; note s50 chr1 YRI
eigonly 0.62); fix results-doc B2 ("GBR/CHB chr22" — s49 ran chr1/2); reconcile contribution 4
and the Fig 6 caption with the new within-population chr1/2 result.

## Verdict

**Almost.** The round-13 experiments are genuine and the core finding (both orderings on real
genomes, with at-chance controls and a replicating in-sample mechanism) is a real contribution.
But the law framing outruns the artifact: a misquoted range hides the chr2 3-pop exception, the
monotonicity claim is contradicted twice, the crossover location is interpolated, and the
size–structure confound is neither disclosed nor controlled. W13-1/W13-2/W13-5 are one-line
text fixes; W13-3's matched rung is a one-parameter rerun (minutes). Land those and this
returns to submission-ready.
