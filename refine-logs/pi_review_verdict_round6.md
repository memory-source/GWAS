# Senior ML Review — Loop 2, Round 1 (2026-10-04)

Reviewer stance: adversarial. Everything below verified directly from artifacts
(paper/main.tex, code/*.py, results/raw/*.json) — executor notes not trusted.
Method: full read of main.tex, common.py + 13 scripts, all 13 raw JSONs; fresh
pdflatex compile; systematic number-match of 66 paper values against JSONs.

## Verification results

- **64/66 checked paper numbers match the raw JSONs** (the 2 "misses" were my
  own rounding tolerance: 0.0921→"+0.09" and 0.1020→"+0.10" are correct).
  Tables 1/2/3, per-ancestor values, CIs, msprime, real-data, s8 functional
  forms, s3b paradox — all match.
- Two real numeric slips: (a) Contribution 2 "+0.003±0.006" vs
  `delta_T2_std = 0.00472` (should be ±0.005); (b) abstract + contribution 1
  "0.85 at full precision" vs F5_full = 0.8446 (mechanism ¶ correctly says 0.845).
- Code implements what the paper claims: 3-way disjoint splits (contiguous
  slices of exchangeable i.i.d. individuals — sound); s7 shuffles within
  super-population before splitting (ancestry confound fixed); maf() is the
  fixed f=mean/2 version (3814/4000 tested ✓); AF-LRT orientation correct;
  PC-distance naive (0.12)/adversarial (0.88) orientations disclosed and both
  supported by artifacts; CV meta-attack is genuine out-of-fold StratifiedKFold;
  controls all real (ref-basis 0.502/0.510, permutation 0.486, T2 refit
  +0.003, T3 attacker-realistic +0.10, eig ratio 1.000, noise 0.496).
- s3b genuinely implemented (computed stopping rules; quality = re-run
  association pipeline; no hardcoded probes). Real 1000G chr22 VCF present
  (205 MB); s7 parse fixed (full split, both alleles).
- Compile: 6 pages, 0 errors, **1 undefined reference** (`sec:discussion` —
  §7 has no `\label`; contribution 4 renders "Section ??"), 3 overfull hboxes.

## Weaknesses (ranked)

1. **Release-unit ambiguity behind "attacker-implementable from the release
   alone"** — the projection-distance attack needs the k×m SNP loading matrix +
   study mean; "PC spectrum" reads as eigenvalues; `release_objects` emits
   eig+scores but NOT loadings, and every attack recomputes loadings from raw
   study genotypes. Claim asserted, not demonstrated against an explicit
   release object. MINIMUM FIX: one sentence naming the released object (PC
   loadings, as federated PCA shares for cross-site projection) + conditional
   scoping; optionally add loadings to release_objects and have T3 consume them.
2. **Operating point of headline experiments unstated** — no n_site/n_ref/n_non
   (300/ancestry → 900/900/900), 4000 SNPs anywhere in the paper. MINIMUM FIX:
   one sentence in §4 + captions.
3. **Uncertainty gaps**: s7 real = single split (seeds=(0,)); boundary and
   dim sweep = 2 seeds, unstated. MINIMUM FIX: state counts; optionally rerun
   s7 with a few split seeds (CPU-hours).
4. **Text/numeric slips**: ±0.006→±0.005; "0.85"→0.84(5); "sign-reversed" is
   wrong (sim 0.12 and real 0.343 naive are the SAME sign — attenuated, not
   reversed); real-data top-200 criterion (|AF−refAF|) differs from simulated
   (association p-values), unstated; Table 2 caption's 5-seed pointer cites
   Table 1 instead of the §6/s2c text. Five text edits.
5. **"Largest intermediate-release channel" comparison class undefined**;
   per-site AF (0.96) is also intermediate; the F/M tier taxonomy that makes
   the claim coherent is never explained yet F-labels appear in Table 1.
   MINIMUM FIX: define taxonomy once / say "largest non-AF diagnostic".
6. **Missingness row near-tautological**: simulated as iid uniform noise.
   MINIMUM FIX: one clause, or generate batch-correlated missingness.
7. **"generally largest at small-to-moderate sites"** — data pattern is clean
   in RATIOS (Δ largest at ratio 3–10, collapses ≥20); swap "sites"→"ratios".
8. **Broken \ref{sec:discussion}** → "Section ??" in PDF. Trivial.
9. Stale docs: EXPERIMENT_RESULTS.md still says real-data PC leak "not
   separately evaluable" (artifact+paper report 0.343); tracker "compiled
   clean" (3 overfull + 1 undefined ref).
10. Code note: auroc() assigns distinct ranks to ties (no average ranks) —
    latent only; all reported scores continuous, no number affected.

## Verdict

Score: **7/10** for a top venue. Every measured number verifies against raw
artifacts; the controls are genuine; the mechanism decomposition (T2/T3) is
the paper's strongest technical asset. What holds it below 8 is positioning
precision (W1/W2), single-split real data (W3), and the small slip list.

**READY: Almost** — submission-plausible after a text-fix pass (W1, W2, W4,
W5, W7, W8 ≈ one afternoon; no re-experimentation required). W3's s7
multi-seed rerun (CPU-hours) would remove the last substantive exposure.
