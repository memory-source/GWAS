# Round-4 Re-Review (autonomous review loop round 2/3) — 2026-10-04

Reviewer: same independent artifact-level review. Read the revised main.tex in full,
re-checked make_figures.py, CLAIMS_FROM_RESULTS.md, results/raw/ (no new artifacts),
recompiled (8 pp, 0 errors, 0 undefined refs), and re-ran the structure measurements.

## Fix verification (against current files, not the author's claims)

- **W1 ADDRESSED (substance)**: both "real spectra are less structured" sentences gone.
  New text (§7 Real data; §8 PC channel) says the smaller real-data magnitude has no
  simple explanation, states the real chr22 panel is MORE structured, and flags the
  shared-super-pop-structure / smaller-in-sample-ness-share cause as a hypothesis.
  I verified the new numbers: **Fst-like ≈0.09 real vs ≈0.05 sim is reproducible**
  (variance of per-pop allele freqs over p̄(1−p̄) with numpy-default ddof=0:
  real 0.0901 vs sim 0.0503; my ddof=1 variant gives 0.135/0.075 — same ratio 1.8),
  and "larger top-10 variance fraction" is true (real 0.262 vs sim 0.144/0.122).
  Direction correct under both variants. Residual: the estimator is unnamed in the
  paper and no script/artifact in the repo backs 0.09/0.05.
- **W2 FIXED**: sentence completed ("holds on real data."); truncation sweep clean.
- **W3 FIXED (substance)**: Discussion now reports the batch variant with the correct
  interpretation (batch independent of membership → coverage-degraded AF, not a new
  channel; membership-correlated composition future work). CLAIMS_FROM_RESULTS.md
  updated ("the current batch-missingness variant is coverage-degraded AF, 0.76, not
  a new channel"). Residual: the parenthetical "(0.76 vs 0.81 at full coverage)"
  mislabels 0.81 — that is the uniform-random-70%-coverage control (my round-1 run,
  0.808/0.822), not full coverage (full = 0.845, artifact-backed via s1c F5_full on
  the same panel config); 0.81 itself has no repo artifact.
- **W4 FIXED**: real-data paired-test sentence matches my round-1 computation exactly
  (+0.043, CI [+0.033,+0.053], 20/20, t p=1.3e-7); derivable from s9_orientation.json
  per-seed data in one line.
- **W5 MOSTLY FIXED**: seed inventory complete (25/20/20/12/20); "sparse-marker" fixed;
  Fig 5(a) now plots ["AF only","PC channel","AF+PC"] with the y-label "(PC orientation
  calibrated)", figure regenerated. **New one-word error introduced**: "peaks near the
  sparse-marker/large-panel corner (at 400 individuals, 500 markers)" — 400 is the
  SMALLEST panel in the S12 grid (400/800/1600); the original "small-panel" was correct.

## New/remaining defects (all text-level)

1. **Line 254: "The calibrated orientation is ``closer'' on 25/25 seeds" is wrong.**
   s10_paired.json: pc_sign = +1.0 on all 25 seeds (sign_agreement 1.0) — the calibrated
   orientation is the members-FARTHER direction, as the paper's own Mechanism paragraph
   says ("members sit systematically farther"). Leftover from the pre-fix version.
   Fix: "farther" (or "agrees with the members-farther direction on 25/25 seeds").
2. **"large-panel" regression** (above). Fix: "sparse-marker/small-panel corner".
3. **"0.81 at full coverage" gloss wrong** (above). Fix: "(0.76, below the 0.845
   full-coverage AF attack on the same panels, Table 1)" — or back 0.81 with a
   dozen-line uniform-coverage control artifact.
4. **Fst-like estimator unnamed / unbacked** (above). Fix: name the estimator
   ("variance of per-population AF over p̄(1−p̄), averaged over SNPs") and add the
   dozen-line script/artifact, or cite it as computed for this comparison.

## Score: **7.5/10** (up from 7.0)

The false explanation is gone and replaced with an honest open hypothesis backed by a
verified-more-structured direction; the real-data increment now has a formal paired
test; the batch caveat is measured and correctly interpreted; the figure is cleaner.
Remaining defects are four one-word/gloss-level text items; no experiments required.

## READY? **Almost**

One text pass (items 1–4, ~15 minutes) from Yes. Item 1 is the one a reviewer with
s10_paired.json would catch (sign error contradicting the artifact); fix it before
submission. Everything else in the paper still verifies against results/raw/*.json.
