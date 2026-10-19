# Round-5 fresh review — Round 4/4 (FINAL)

Reviewer: pi (same method). make_figures.py fig7a, regenerated fig7 PDF, and
main.tex verified; residual-claim grep; fresh compile (12 pp, 0 errors,
0 undefined refs).

## Verification of the four claimed fixes — all applied

1. **Fig 7a band fix: APPLIED and now correct.** axvspan(0.09, 0.16) on the
   nominal axis; the in-figure label reads "measured 0.06 (crossover) → 0.10
   (real 1000G)" — it explicitly says "measured", so the band is correctly
   presented as the nominal-F_ST equivalents of the measured endpoints
   (0.06/0.63 ≈ 0.095; 0.10/0.63 ≈ 0.16, per s21's fst_like ratios 0.62–0.68).
   The code comment states measured ≈ 0.63 × nominal. fig7_governance.pdf
   regenerated after the edit (mtime 03:20 > 03:19). The band now contains both
   labeled points' nominal equivalents — figure, label, comment, and axis agree.
2. **Fig 7a caption: APPLIED.** "crosses zero at a measured structure ≈0.04–0.06;
   the shaded band spans the nominal-F_ST equivalents of that crossover and of
   the (higher) real-1000G structure level (measured ≈0.10, nominal ≈0.16), so
   the real panels sit beyond the crossover" — both units stated, arithmetic
   correct, consistent with §7 and the s21 JSON.
3. **§7 merge: APPLIED.** The "best read in measured terms" line is folded into
   the crossover sentence ("at nominal F_ST=0.08 the realized between-group
   statistic is ≈0.05, i.e. measured ≈ 0.63× nominal"); the paragraph is coherent
   and no longer redundant.
4. **Crossover range: APPLIED in the abstract** ("crosses at measured structure
   ≈0.04–0.06, below the real-data level") — matches the s21 data (measured
   crossover ≈0.062 at n800, ≈0.043 at n1600). Contribution 5 still says
   "crossing at measured ≈0.06" (the n800 value) — factually correct but one
   step behind the abstract/caption in precision (see cosmetic note).

Residual-claim sweep: zero matches for all previously flagged phrases. Regression
check: headline numbers, bibliography, governance numbers, canonical scoping,
IBS descriptions — all unchanged and hold. Compile clean (12 pp).

## Score: 9/10

## Verdict: YES — submission-ready

No remaining defects found. Every checked number matches the raw JSONs; all
interpretive claims now match what the code measures; scoping is honest
throughout (PC dominance simulation-scoped, real-data attenuation reported,
prototype honestly not-MPC, sole-egress binding argued not built).

## Remaining weaknesses + minimum fixes

1. **Cosmetic (non-blocking, one word):** contribution 5 says "crossing at
   measured ≈0.06" while the abstract and Fig 7a caption say "≈0.04–0.06".
   Unify to "≈0.04–0.06" (or "≈0.06 at n=800, ≈0.04 at n=1600") for consistency.
2. **Cosmetic (non-blocking):** the code comment says the endpoints "map to
   nominal ~0.10 and ~0.16" while the band is drawn at 0.09–0.16 (the 0.06
   point maps to ≈0.095) — a rounding-level comment/code mismatch inside
   make_figures.py only; invisible in the paper.

## Assessment

Across the four rounds: 8 → 8.5 → 9 → 9 with all blockers resolved. The
governance-layer round closed all five prior limitations with real experiments
(Fst sweep, admixture stress test, working manifest/egress/ledger prototype with
a real violation matrix and composition measurement, the paradox as a controlled
falsifiable result with a null control and negative side, and the IBS bound with
a planted control). The three interpretive defects found in round 1 (inverted
Dirichlet reading, nominal/measured conflation, "phased" mislabel) were fixed
correctly by rounds 2–3, and the final figure-axis detail in round 4. The paper
is submission-ready as a measurement paper.
