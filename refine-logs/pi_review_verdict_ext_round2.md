# Review verdict — extension round 2 (fix verification after ext-round-1 8.5/10 READY)

**Reviewer**: pi (independent; all changed files read/verified directly)
**Scope**: W1–W4 fixes + refreshed s53 (9 real cells) + paper/docs propagation
**Verdict: READY — 9/10**

## Fix verification (all confirmed in the files)

- **W1 (tracker wording) — FIXED.** Tracker L1 now reads "ladder falls
  monotonically and crosses on chr22 (+0.120/+0.042/+0.007/−0.026) and chr1
  (+0.058/+0.002/−0.010/−0.018); chr2 stays positive, near-parity at the top
  rung (+0.109/+0.059/+0.018/+0.006)". Paper: "crossing locations on chr22/chr1
  (n_study≈260–400)". CLAIMS C6 and EXPERIMENT_RESULTS L1 scoped identically.
- **W2 (seeds paragraph) — FIXED.** Parenthetical now includes "the fixed-basis
  probe / second-haplotype / real size-ladder studies 5 splits or seeds per
  cell" (main.tex ~line 263).
- **W3 (per-rung noise) — FIXED and accurate.** "The top rungs sit within split
  noise of zero; chr2 stays positive" — verified against the JSONs: chr22
  per120 +0.007±0.026, chr1 per80 +0.002±0.059 / per120 −0.010±0.022, chr2
  per160 +0.006±0.023 (all within 1 split-sd of zero).
- **W4 (s53 coverage) — FIXED, and extended cleanly.** Script now runs
  chroms=(22,1,2) × sups=(EUR,EAS,AFR) → 9 real cells; JSON refreshed 03:24
  with script mtime 03:22 (edit preceded the run — round-1's post-run-edit
  pattern is gone). Decisive integrity check: the sim cell AND the four
  carried-over real cells are value-identical to the round-1 JSON (same seeds,
  deterministic, additive-only extension). New cells agree with the s51 pooled
  arm on the same pools within ~0.02 (different calibration masks explain the
  small deltas). Fixed bases 0.484–0.530 (worst deviation 1.5 split-sd,
  chr1:AFR dct 0.526±0.017 — unremarkable across 45 fixed-basis cells and
  disclosed inside the quoted range); full loadings 0.787–0.903.
- **Propagation — COMPLETE.** Paper quotes 0.46–0.53 + chr22/chr1/chr2×EUR/EAS/
  AFR in all three s53 spots (abstract, mechanism, discussion); TRACKER L4
  (".46-.53 … .79-.90"), RESULTS L4, CLAIMS C16 all updated; no stale
  "0.48–0.52"/"0.81–0.89" anywhere; s51/s52 passages untouched.
- **Compile — independently verified**: pdflatex ×2, 0 errors, 0 undefined
  references, 21 pages (matches the claim).

## Remaining weaknesses (both trivial, neither blocking)

- **N1 (low, one-number attribution):** the quoted "0.46–0.53" range for
  "fixed, data-independent bases" (abstract + mechanism paragraph) is the union
  with the random basis: fixed bases alone span 0.484–0.530; the 0.464 endpoint
  is the random basis (chr2:EAS rand 0.464). The discussion's wording ("an
  eigenvalue-only (random-basis) projection, and fixed data-independent bases …
  are at chance (0.46–0.53)") is exactly correct; the other two spots attribute
  the union range to fixed bases alone. Conclusion unaffected (all 45
  fixed-basis cells in [0.484, 0.530], full-loadings contrast 0.79–0.90).
  MINIMUM FIX: in the abstract and the mechanism paragraph, either quote
  "0.48–0.53" for the fixed bases or extend the subject to "random-basis and
  fixed bases" to match the discussion. Two one-line edits.
- **N2 (cosmetic):** s53_eigbasis.py docstring still says "(EUR/EAS on chr22
  and chr1)" while the code/JSON now cover 3 chromosomes × 3 super-pops
  (same docstring-lag nit class as round 1; JSON self-documents via
  chroms/sups). MINIMUM FIX: one-line docstring edit.

## What plainly holds

The extension round did exactly what it claimed, the fixes are applied
consistently everywhere, the refreshed s53 is a clean additive re-run whose
unchanged cells reproduce bit-for-bit, and the compile is clean. The paper's
disclosed structural limits (1000G-only real data; small real n for the
within-population PC>AF cells; semi-synthetic large-n arm; 5-split ladder
uncertainty, now explicitly quantified as "within split noise") remain as
honestly scoped design properties, not defects.

**Score: 9/10 — READY.** (Round 1: 8.5. The two residual items are one-line
cosmetics with zero claim impact; withholding the last half-point for the
paper's inherent, disclosed scope limits rather than any remaining defect.)
