# Round-15 Review, Loop 2, Round 3 (re-review after R15-2's 7.5/10) — run_20261010_a1f2c3d4

**Reviewer**: pi. **Date**: 2026-10-10.
**Scope**: verification of the five fixes claimed against R15-2's W1–W5; regression check of the
edited paragraphs; compile state.

## Fix verification (all checked in the changed files, not trusted)

- **W1 (s53 into the body) → FIXED (with one new nit).** Mechanism paragraph (main.tex:415–419)
  now reads "…the random-basis result is a specific (eigenvalue-only projection) ablation, not a
  proof that no eigenvalue-dependent attack exists, and fixed, data-independent bases---including
  one weighted by the released eigenvalue spectrum---are also at chance (0.48–0.52 …), so a purely
  eigenvalue-informed fixed basis does not recover the leak." Discussion (:1172–1175) now reads
  "an eigenvalue-only (random-basis) projection, and fixed data-independent bases including an
  eigenvalue-weighted one, are at chance (0.48–0.52) while the full study-fitted loadings recover
  the leak." Both sentences are accurate against `s53_eigbasis.json` (0.484–0.517). The abstract
  (:73–75) is unchanged and consistent.
  **New nit N1**: the added mechanism cross-reference `Section~\ref{sec:replication}` resolves, but
  sec:replication contains no s53/fixed-basis content (0.48–0.52 appears only at :75, :418,
  :1174); a reader following the pointer lands on the *random-basis* ladder numbers (0.51–0.516).
  Fix: drop the pointer or retarget (the mechanism sentence itself is now the body report).
- **W2 (stale docs) → FIXED.** `EXPERIMENT_RESULTS.md` L1 bullet now carries the joint-size-effect
  wording, the interpolation label, and the "zero-structure real cohort at that size is not
  available inside 1000G" admission; L5 bullet now "the full 48000-row fit was not run" (no GB
  claim); tracker L4 row ".48-.52" and "cosine". All three match the paper.
- **W3 (claims doc) → FIXED (with one new nit).** Evidence line includes all five round-15 JSONs;
  the C6 prose bullet now reports the fixed-composition ladder with numbers that match
  `s51_real_lowstruct_chr{1,22}.json` exactly (+0.120/+0.042/+0.007/−0.026;
  +0.058/+0.002/−0.010/−0.018; ~260–400 interpolated; pooled arm 167–220 PC>AF); the
  "do not support" entry now records the cap explanation and "a joint-size crossing, not a
  zero-structure one"; new table rows C14/C15 are accurate in content.
  **New nit N2**: the new table row "C14 — fixed, data-independent bases…" collides with the
  existing prose bullet "C14 (round 23): as-released ordering…" — two different claims now share
  the ID C14 in one document. Renumber the new rows C16/C17 (or renumber the prose C14).
- **W4 → FIXED.** "Two further leave-self-out haplotype statistics" (:1053); rest of the sentence
  and all eight s52 numbers unchanged and correct.
- **W5 → FIXED.** s38 comment now "af_full has the higher AUROC but the LOWER eps" — correct
  (matches the JSON inversion note).
- **W6 → open, as declared.** 7 overfull boxes (largest 88.7pt), cosmetic.

## Regression check

The mechanism and Discussion edits are additive; all previously verified numbers in the edited
paragraphs are intact (k-ladder 0.53/0.65/0.74/0.88; rand-basis 0.51; conditioned-AF paragraph;
s52 numbers). Compile: 21 pages, 0 errors, 0 undefined references; PDF (02:07) newer than source
(01:58).

## 1. Score: 8 / 10

## 2. Verdict: **READY**

The manuscript now states exactly what its artifacts show at every position checked across
R15-1/R15-2/R15-3, including the body positions that argue its punchline. The residual defects are
copyedit-class and do not touch any claim, number, or experiment. The honest ceiling is set by the
disclosed scope bounds (zero-structure large-n cohorts unobtainable inside 1000G — now with a
measured cap explanation; semi-synthetic scale arm; tested-attack scoping), not by defects.

## 3. Remaining weaknesses (non-blocking; minimum fixes)

- **N1 (copyedit).** main.tex:418 cross-reference `Section~\ref{sec:replication}` points to a
  section that does not discuss the fixed-basis result. Fix: delete the pointer or change it to
  the Discussion reference.
- **N2 (copyedit).** CLAIMS_FROM_RESULTS.md: table row C14 collides with prose bullet C14.
  Fix: renumber the two new table rows C16/C17.
- **W6 (cosmetic, open).** Two overfull hboxes at 76–89pt. Fix: rewording or \linebreak in the s52
  paragraph and Limitations; optional.

## 4. Minimum fixes before submission

N1 + N2 (a five-minute pass). W6 optional. No further experiments, no re-verification round needed.
