# Round-2 Adversarial Review — Selection-Input Spectral Tomography (deep loop, round 2)

**Reviewer**: pi (senior-ML-reviewer protocol, adversarial; artifacts read and re-run directly, nothing trusted from the executor)
**Date**: 2026-10-16
**Round-1 verdict**: 6.5/10, ALMOST (W1–W7 listed in `pi_review_verdict_deep_round1.md`).

## Score: **7/10 as-is. Verdict: ALMOST** (ready after ~half a day of text corrections and two factual fixes; **no new experiments required**; one optional 15-minute cell).

---

## A. What I verified myself this round

1. **Exact reproduction under the *current* code — the new code-postdates-results concern is resolved.**
   The s61 code was edited at 18:00 (adding the `qs` parameter) *after* the s61/s62 JSONs
   (17:22). I re-ran `s61.witness` + `s61.cell` (within:EUR, chr22, seed0) with the current
   code: every numeric leaf matches the JSON to ≤1.4e-8 (BLAS noise; the largest diffs are in
   degenerate `target_err_inside` terms), the block-swap sequence is **bit-identical**
   ([0,0,22,28,4,12,20,26,28,12,24,22]), and all TPR/no-release/precision-probe values match.
   The 18:00 edit is behavior-preserving on the default path.
2. **W1 (TPR@FPR) — genuinely fixed at code level and in the chr22 tables.** `metrics()`
   returns `membership_tpr_fpr01/05`; `summarize()` aggregates them plus
   `baseline_tpr_fpr05_single_release` and `spec_cond_tpr05_mean`; s63's
   `_membership_auc_from_operator` returns `(auc, tpr05, sign)` and its per-seed JSON carries
   TPR values. Verified against `s61_spectral_tomography_chr22.json`: random q8 TPR@5% =
   1.000/0.580/1.000 (EUR/3pop/EAS), single release 0.675/0.139/0.616, conditioned single
   release 0.916/0.299/0.899 — exactly as the results md now reports. The chr22 gap tables
   all carry the new column with correct values.
3. **W2 (saturation) — new cell run and verified; claim needs scope before integration.**
   `code/s64_saturation.py` faithfully reuses `s61.cell` (added `qs` kwarg; the
   precision/sci-cost probes are correctly guarded to the default `qs`). The JSON matches the
   executor's numbers exactly: random q8→q16 gap +0.0151→+0.0153 (CI [0.0088,0.0219], 5/5
   positive), visible share 0.853→0.927, rank 72→130, recon err 0.521→0.376; adaptive
   +0.0152→+0.0153; passive flat. **But**: (a) s64 appears in *neither* the results md (mtime
   17:55 < s64 18:15) nor the claims file; (b) the saturation is demonstrated only on the
   cohort where the attack is already ceiling-censored — EUR random q8 has TPR@5%=1.000 and
   AUC 1.000, so a flat q8→q16 gap is partly "cannot exceed perfect". The cohort with
   headroom (3pop: AUC 0.907, TPR@5% 0.58 at q8) was not run at q16. "The q-effect on
   membership is saturated, measured rather than extrapolated" is an honest *EUR-cell*
   statement and an overclaim as a general one.
4. **W3 (no-release control) — implemented exactly as specified, values verified, but one
   reported number is an artifact (see F1).** The control (cal-member PCA spans at d=k and
   d=qk, plus a random k-dim Gaussian span, scored through the same orientation path) is in
   `cell()`, averaged in `summarize()`, present in **all four** s61/s62 JSONs, and has its own
   results-md section. chr22 means verify: EUR 0.466/0.457/0.490; 3pop 0.500/0.494/0.485;
   EAS 0.490/0.450/**0.579**.
5. **W4 (adaptive reframe) — half done.** The "designed measurement" overclaim is gone from
   claims+md; D6 is stated as a negative result in both. **The claimed stronger-threat
   sentence ("random AF-preserving swaps suffice; a defender cannot rely on detecting
   targeting") appears in neither file** — grep over both documents returns nothing.
6. **W5 (heterogeneity reframe) — done in the md.** Support item 4 foregrounds 3pop +0.123,
   states the within-population increments are "small", and the n-sensitivity section is
   retained with correct numbers (I re-verified every s62 cell: +0.011/+0.015 at n120,
   +0.034/+0.051 at n180, CIs match, honest120 ⊂ honest180 with shared 120/80/60/40 roles).
   The fragile EAS-q4 +0.004 rows are retained but with honest CIs (chr1's CI includes zero,
   reported as 4/5) — acceptable, see R6.
7. **W6 (interface scoping) — NOT done.** The scoping sentence (claim covers interfaces
   returning joint top-k loadings+eigenvalues per re-fit; SF-PCA-style selective decryption)
   exists only in the *pre-existing* direction document (`DEEP_DIRECTIONS_ZH.md`), not in the
   claims file or results md.
8. **W7 (nits) — partially done.** `Zh` param removed from `metrics()` ✓ (verified). The
   histogram/AF/variance statement is explicit, but HWE/missing preservation is still never
   spelled out; the "span-restricted (C_H P_Q), never full C_H" phrasing was not added;
   `_boot_ci` still hard-codes seed 7 undocumented (fixed bootstrap seed is fine practice —
   but the promised note is missing).
9. **s63 reconciliation** — every md range matches the JSON (af 0.615–0.775, af_cond
   0.855–0.924, single 0.879–0.951, PLINK-conditioned 0.976–0.994, stacked q8 0.994–1.000;
   gap +0.0140 [0.0074,0.0206], 5/5; identity residual 6.2e-6 over 40 exports; marker-order
   check true).

## B. New findings this round

**F1 (the important one) — the EAS random-span no-release control (0.579) is a preprocessing
artifact, and the md sentence "All at or below chance" is factually wrong.**
I probed the EAS random_k = 0.579 (5/5 seeds above 0.5: 0.534–0.634) with 20 *fresh* random
spans per seed: the oriented control centers at 0.53–0.60 on EAS but 0.508 on EUR, and the
*unoriented* AUROC is ≈0.41 — honest members systematically score **lower** than nonmembers
on *any* fixed span. Decomposition: out-of-`mu` rows have 5× larger z-norms on EAS
(nonmem RMS z≈4.7 vs honest 0.94; the effect flips exactly when the `mu` panel changes).
Mechanism, confirmed in code: `within:EUR` is loaded via `s48.load_for(("EUR",),…)` (MAF
filter on EUR), but `super:3pop`/`within:EAS` load via `s51._sup_label(SUPERPOP,…)` →
`_parse_vcf(..., sups=all)` → **MAF filter on all 2504 samples**. The EAS panel therefore
contains globally-common, EAS-rare markers: on chr22-EAS seed2, **392/3271 markers have
`mu==0` exactly** over the 240 honest+partner rows (113 more with mu<0.01), the variance-
standardization scale hits its 1e-6 clip, and any out-of-`mu` carrier gets |z| up to 2000
(13/80 nonmember rows have max|z|>20). A random Gaussian span picks these up; the
orientation step flips the resulting anti-membership signal into +0.08 above chance.
**Decisive control**: re-running the full EAS cell with |z|≤5 clipping moves
`random_k` 0.579→**0.463** (chance, where it should be) while leaving every headline number
intact — af_cond 0.915→0.915, spec_cond 0.987→0.986, **q8 gap +0.0126→+0.0140, all 5 seeds
positive** (unpatched rerun reproduces the published JSON exactly). Blast radius: this one
control cell and the md sentence containing it. chr1 controls (0.477–0.519), chr22
EUR/3pop controls, all cal-PCA controls (spans from in-`mu` rows have exactly zero weight
at mu==0 markers), s62 (EUR panel), s63 (EUR panel) and s64 (EUR panel) are clean. The
leakage claims are *unaffected* — if anything the corrected control (0.463) strengthens
"the release-derived span carries the signal".

**F2 — "adaptive−random CI includes zero at every q in EUR/3pop" is false as written.**
chr1 3pop q4: +0.0097, CI [+0.0031,+0.0185] — excludes zero (it vanishes by q8: +0.0043,
CI includes 0). chr22 EUR q8: [+0.00000,+0.00030] — degenerate (4 seeds tied at ceiling).
The substantive D6 conclusion (no *material* adaptive advantage anywhere; max +0.010) is
right; the sentence needs one clause, not a reanalysis.

**F3 — integration debt.** `EXPERIMENT_AUDIT_DEEP.{md,json}` (05:26) predates every round-2
artifact; the claims file (09:44) predates the afternoon re-runs (its D1–D7 numbers still
reconcile — I checked — but it contains none of the round-2 additions); s64 is absent from
both paper-facing documents. My own replay (A1, A3) closes the integrity gap for this round.

## C. What is correct and should be said plainly

- Reproduction is exact and the 18:00 edit is harmless; no honest-matrix access in the
  attack path (re-inspected `adaptive_input`, `metrics`, `_membership_auc_from_operator`,
  the no-release block); negative results are still reported, not buried (q1 losses, D6,
  D7, the 4-dp boundary).
- The TPR@FPR addition (W1) genuinely improves the paper: "all members recovered at 5% FPR
  vs ≤0.92 for the strongest single-release attacker" is the right headline and is now
  measured, not asserted.
- The no-release control's *cal-PCA arms* — the objection that actually matters ("isn't this
  just knowing some members?") — are sound and clean in every cohort (0.41–0.52).
- The s64 machinery is honest (same seeds, same code path, no new degrees of freedom), and
  its coverage-side finding (share 0.93, recon err 0.38 at q16 while membership is flat) is
  a real, useful measurement.

## D. Remaining weaknesses, ranked, with minimum fixes

| # | Weakness | Minimum fix |
|---|----------|-------------|
| R1 | F1: EAS `random_k`=0.579 reported as "at or below chance"; control contaminated by global-MAF/scale-clip artifact | Correct the sentence; re-emit the three no-release controls with \|z\|≤5 (or drop markers with cohort-mu==0) and report the corrected values (EUR/3pop/EAS ≈ 0.49/0.48/0.46); state explicitly that all headline numbers are unchanged (my rerun: gap +0.0140 clipped vs +0.0126 raw) |
| R2 | W4's threat sentence and W6's scoping sentence were claimed but never written | Add one sentence each to claims + results md: "uniformly random AF-preserving swaps suffice — targeting is unnecessary at q≤8" (appendix keeps D6); "the claim covers interfaces returning joint top-k loadings and eigenvalues per re-fit (PLINK 2.0 `--pca allele-wts`; SF-PCA-style selective decryption), not all federated GWAS pipelines" |
| R3 | F2: adaptive−random CI wording contradicted by chr1 3pop q4 | One clause: "never exceeds ~0.01; the only CI excluding zero (chr1 3pop q4, +0.010 [0.003,0.018]) vanishes by q8" |
| R4 | s64 absent from the paper; saturation claim would be an overclaim if integrated as-is | Integrate with scope: "EUR chr22 only, where TPR@5% is already 1.000 at q8; coverage keeps rising (share 0.93) while the membership gap is flat q8→q16". Optional (recommended, ~15 min): the same cell for 3pop, the one cohort with headroom |
| R5 | TPR columns missing from chr1/s62/s63 md tables; claims file lacks the practical reading; HWE/missing, span-restricted phrasing, `_boot_ci` seed note still unwritten; audit stale | Text edits + refresh `EXPERIMENT_AUDIT_DEEP` against the 17:22–18:15 JSONs |
| R6 | EAS-q4 +0.004 (chr22, 5 seeds, CI [+0.001,+0.006]) still carries the "appears from q≈4" generalization | 15–20 seeds on within-population arms (~37 s/cohort-seed) or drop EAS from the q≈4 sentence (3pop +0.047 and EUR +0.011 carry it) |

## E. Verdict

**ALMOST — closer than Round 1.** The three substantive round-1 asks (TPR@FPR, no-release
controls, a saturation cell) are genuinely implemented and their numbers verify against the
raw JSONs; my attempted break of the new machinery failed everywhere except one control
cell, and that failure (F1) turned out to strengthen, not weaken, the leakage claim. What
remains is text: two sentences that were claimed but never written (W4/W6), two factual
corrections (F1's "all at or below chance", F2's CI clause), and the s64 integration with an
honest censoring scope. With R1–R4 applied this is a credible bounded measurement paper
(7.5/10) for a security/measurement venue; it is still not a methods contribution and should
not be framed as one.

---

## Memory Update (for next round)

- **Round-2 state**: score 7/10, ALMOST. W1/W2/W3 substance verified genuine; W4 threat
  sentence, W6 scoping sentence, W7 phrasings claimed-but-absent; s64 not yet integrated.
- **F1 artifact (key fact)**: EAS/3pop panels are global-MAF-filtered (`s51._sup_label` →
  `_parse_vcf(sups=all)`) while EUR is EUR-filtered (`s48.load_for(("EUR",))`); chr22-EAS has
  392/3271 markers with cohort-mu==0 → scale-clip 1e-6 → |z|≤2000 for out-of-mu carriers →
  no-release `random_k` = 0.579 (md wrongly says "at or below chance"). Clipped-|z|≤5 rerun:
  random_k 0.463, q8 gap +0.0140 (survives; 5/5), spec_cond/af_cond unchanged. chr1, s62,
  s63, s64 all clean (EUR-filtered or mixed panels). Reviewer scripts:
  `scratch/_rv2_repro_and_probe.py`, `scratch/_rv2_eas_decomp.py`,
  `scratch/_rv2_eas_clip_rerun.py`.
- **Reproduction**: current s61 code (18:00, `qs` param) reproduces s61 chr22 EUR seed0 to
  ≤1.4e-8, bit-identical swaps — code-postdates-results resolved for this round.
- **Watch next round**: (1) the corrected no-release table; (2) the two missing sentences;
  (3) s64 integrated with the EUR-only/ceiling-censoring scope (or a 3pop q16 cell);
  (4) adaptive−random CI clause (chr1 3pop q4 +0.010 [0.003,0.018]); (5) audit/claims
  refresh. Do not accept "saturated" as a general claim from the EUR cell alone.
