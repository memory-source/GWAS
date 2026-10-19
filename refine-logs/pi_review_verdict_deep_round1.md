# Round-1 Adversarial Review — Selection-Input Spectral Tomography (deep loop)

**Reviewer**: pi (senior-ML-reviewer protocol, adversarial; artifacts read directly, nothing trusted from the executor)
**Date**: 2026-10-16
**Artifacts verified**: `DEEP_DIRECTIONS_ZH.md`, `EXPERIMENT_PLAN_DEEP_20261016.md`,
`CLAIMS_FROM_RESULTS_DEEP_20261016.md`, `code/s61_spectral_tomography.py`,
`code/s62_sample_size.py`, `code/s63_plink_export.py`, `code/common.py`,
`code/real_multi.py`, `code/s48_real_structure.py`, `code/s51_real_lowstruct.py`,
all four raw JSONs, `EXPERIMENT_RESULTS_DEEP_20261016.md`, `EXPERIMENT_AUDIT_DEEP.md/.json`,
`review-stage/DEEP_REVIEW_20261016.md`.

## Score: 6.5/10 as-is → 7/10 after the fixes below. Verdict: ALMOST (ready after ~1 day of cheap additions; no new experiments required beyond one optional saturation cell).

---

## A. What I verified myself (not taken from the audit or the codex thread)

1. **Exact reproduction.** Re-ran `s61.witness` + `s61.cell` (within:EUR, seed0, chr22)
   with the current code: all membership/coverage/baseline metrics match the JSON to
   ≤ 4.5e-9; the block-swap sequences are bit-identical
   (`genotypes_changed_per_release` [0,0,22,28,4,12,20,26,28,12,24,22] reproduced).
   The 08:25 code mtime postdating the 08:13 JSON is therefore harmless — current code
   reproduces current results.
2. **Full number reconciliation.** Every table in `EXPERIMENT_RESULTS_DEEP_20261016.md`
   and every row of the claim map reconciles with the raw JSONs: s61 chr22 witness
   (max residual 3.820e-11 over 15 entries, AF diff exactly 0.0, hist identical,
   het shift ≤ 4.9e-3), all baseline/gap/CI tables, chr1 replication, s62 nested
   ladder (role sharing verified: partner/nonmember/ref/shadow identical, honest120 ⊂
   honest180), s63 (40 exports, max relative residual 6.242e-6, order check true,
   +0.0140 [+0.0074,+0.0206], 5/5). No discrepancy found anywhere.
3. **No information leakage into the attack.** Code inspection: `adaptive_input`
   reads only partner data, shadow panel, past released V's, and the public
   center/scale (which is a deterministic function of the released pooled AF, so it
   grants nothing beyond the release). `metrics` uses `C_H` only for evaluator-side
   coverage denominators, never in the attack score. The membership score is a plain
   projection norm; orientation is fixed on the calibration half, evaluated on the
   disjoint held-out half.
4. **The identity is load-bearing — my adversarial control.** I ran the attack with
   the tomography identity *removed* (stacking the raw released eigenvectors V_t and
   scoring projection onto their union, 5 seeds × 3 cohorts, chr22). Raw-V q8 gaps vs
   the strongest passive baseline are **negative in all three cohorts**
   (EUR −0.004 [4/5], 3pop −0.012 [2/5], EAS −0.050 [0/5]) while the identity-based
   stack is positive (+0.015/+0.123/+0.013, 5/5 each). The mapping
   Y_t = V_tΛ_t − K_tV_t is therefore not decoration; the capability genuinely comes
   from measuring through the honest scatter, not from "PCA outputs have big spans".
5. **Missing no-release control, run by me.** An attacker holding the 60 calibration
   members' full genotypes but *no release*, scoring projection onto the cal-PCA span
   (d=10/40/80), gets AUC 0.47/0.42/0.43 — at or below chance (ref-span and
   random-span controls: 0.59/0.41). So the release-derived span is qualitatively
   different from anything the cal grant provides: the membership signal requires the
   targets' own directions to have entered the released object. This *strengthens*
   the leakage claim and must be in the paper (see W3).
6. **Ceiling check.** With the paper's orientation convention, per-seed spec_cond is
   0.976–0.994 (matches JSON; my first quick script's seed-3 collapse was my own bug,
   not the paper's). The q8 attack hits TPR@FPR=0.05 of **1.000 on 5/5 seeds** vs
   0.83–0.98 for the strongest passive baseline: the practical effect is larger than
   the AUC delta suggests, and AUC at 1.000 is censored (see W1).
7. **Timestamps coherent**; audit (05:26) predates the final s61 JSONs (07:07/08:13),
   but my own replay closes that gap.

## B. Ranked weaknesses and minimum fixes

**W1 (most important) — Effects are reported only as AUC, and AUC is ceiling-censored.**
Within-population headline increments (+0.015 EUR, +0.013 EAS) sit on a saturated
scale (attack AUC 0.999–1.000, baseline 0.985). A reviewer cannot tell whether +0.015
is meaningful; TPR@FPR (already implemented in `common.tpr_at_fpr`, never used) shows
it is: 1.000 vs 0.83–0.98 at FPR 0.05. The paper currently *undersells* its own effect
and invites the "so what" rejection.
*Minimum fix*: add TPR@FPR ∈ {0.01, 0.05} columns to every gap table (s61, s62, s63,
chr1) and state the practical reading ("all members recovered at 5% FPR vs ≤98% for
the strongest single-release attacker"). Half a day, no new runs.

**W2 — The core capability question is answered at a single point of the design space.**
One k=10, one block length 20, one candidate budget, q ≤ 8, chromosome prefixes. At
q8 only 66–73 of the 120 honest dimensions are visible (scatter err 0.29–0.53): the
q-scaling is nowhere near saturation, so "how many independent measurements can it
design" — the paper's own stated question — is extrapolated from a monotone prefix.
*Minimum fix*: either one added saturation cell (q=16, or k∈{5,20}, one cohort,
~1–2 h) or an explicit scoping sentence in abstract/claims limiting the capability
claim to the tested regime (q ≤ 8, k = 10). The text fix is acceptable; the cell is
better.

**W3 — The baseline ladder lacks the no-release control.**
The obvious reviewer attack is "isn't this just knowing some members?" My control
(W-A5 above) shows the answer is no (cal-genotype-PCA spans are at/below chance), but
the paper never runs it, so the absolute level of spec_cond (0.985) is uninterpretable
in-text.
*Minimum fix*: add the cal-genotype-PCA-span row (d = k and d = qk) plus the
ref/random-span negative controls to `af_baselines`; ~15 lines, reuses existing
roles. Report "release-derived span ≫ known-member-genotype span".

**W4 — The adaptive result is a dead planned contribution that is still framed as machinery.**
D6 (adaptive ≈ random) leaves the direction document's "designed measurement" story
with no operational targeting: target residual only 0.99 → 0.83 at q8, and
adaptive−random CI includes zero everywhere that matters. Keeping adaptive in the
headline framing invites "your only novel algorithm does nothing".
*Minimum fix*: reframe explicitly as a *stronger threat statement* — "uniformly
random AF-preserving swaps suffice; a defender cannot rely on detecting targeted
behavior" — and demote adaptive to an appendix negative result (the results md
already half does this; the claims/abstract do not).

**W5 — Effect-size heterogeneity is under-narrated.**
The within-population increments are +0.013–0.015 (and EAS q4 +0.004 with CI
[+0.0014,+0.0061] on 5 seeds — fragile), while 3pop is +0.123 and n=180 is +0.051.
The paper's one-sentence summary ("small increment") buries the one result with
practical bite and lists one that will not survive scrutiny.
*Minimum fix*: 15–20 seeds for the within-population arms (≈37 s per cohort-seed) or
drop the sub-0.005 EAS-q4 claim; foreground 3pop and the n-sensitivity as the
quantitative headline.

**W6 — Interface generality is asserted from one software package.**
PLINK 2.0 only; no federated pipeline; roles drawn from one 1000G release. Already
acknowledged, but the abstract-level wording ("final outputs as measurement oracles")
is broader than the evidence.
*Minimum fix*: one scoping sentence: the claim covers interfaces that return the
joint top-k loadings and eigenvalues per re-fit (e.g., SF-PCA-style selective
decryption), not all federated GWAS pipelines.

**W7 — Minor code/report nits.** `metrics()` takes an unused `Zh`;
`_boot_ci` hard-codes seed 7; witness HWE/missing preservation is implied by exact
histogram equality but never stated as such in the report; passive q2–q8 rows are
definitionally identical repeats (documented, but say "definitionally" once); the
scatter-recovery claim should always carry "span-restricted (C_H P_Q), never full
C_H" next to the numbers (the md says "no honest-cohort reconstruction" once, in the
negative-results list).

## C. What is correct and should be said plainly

- The identity (A) is one-line algebra, and the paper does not oversell it; the
  non-trivial content — survival through a real PLINK 2.0 export at printed
  precision (6.2e-6 relative over 40 exports, marker order checked per marker) — is
  genuine and I verified the read-back path (ALT-row interleave is order-validated).
- The AF-preservation lemma (B) is exact: AF change 0.0, histogram equality on every
  release, het shift ≤ 4.9e-3 across 15 witness cells on two chromosomes.
- The negative results are reported, not buried: q1 loses everywhere (−0.076/−0.129/
  −0.072), D6 adaptive = no, D7 export precision = no for scatter recovery. No
  selective omission found anywhere in this round (contrast with the previous loop's
  round-3 finding).
- No fabrication, no post-hoc code/result mismatch, no honest-matrix access in the
  attack path.

## D. Verdict

**ALMOST.** The measurement core is sound and honestly reported; I could not break
it after trying the four attacks a hostile reviewer would try (no-release
controls, identity-ablation, orientation exploitation, ceiling analysis). It is not
yet submission-ready *as presented* because W1 (AUC censoring) and W3 (missing
no-release control) leave the two most predictable reviewer objections unanswered,
and both fixes are under a day of work with zero new data collection. With W1, W3,
W4 and the W2 scoping sentence applied, this is a credible bounded measurement paper
(7/10) for a security/measurement venue; it is not a NeurIPS/ICML methods
contribution and should not be framed as one.
