# Cross-model re-review verdict — round 2 (2026-10-02)

Reviewer scope: hostile-but-fair. I verified every revised number against
`results/raw/*.json` (all match the bundle and `EXPERIMENT_RESULTS.md`), read
the new scripts (`s1b_disjoint.py`, `s2b_intermediate.py`, `s5_boundary.py`,
revised `common.py`), and — because the Δ+0.108 question could not be settled
by reading — **ran a diagnostic** (`results/raw/_s2b_diagnostic.py`) that
reproduces S2b exactly (base+M2 per-seed 0.850/0.861/0.871) and adds four
controls. All verdicts below are grounded in code and experiment, not prose.

## Fix ledger (claimed vs verified)

| claimed fix | verified? | evidence |
|---|---|---|
| Disjoint non-members | **YES** | `split_population3` three-way split; s1b/s2b/s5 score `nonmem`, ref used only for `r_af` |
| Bootstrap CIs | **YES** (S1b TPR@1% only) | `bootstrap_ci` in s1b; none for S2b despite its docstring promising "permutation decomposition and bootstrap CIs" |
| S3 redesigned (non-degenerate search, measured quality) | **NO** | `s3_paradox.py` (mtime 03:39) and `s3_paradox.json` (03:52) are byte-identical to round 1: `probe_sequence()` still returns hardcoded `[1,2,3]`/`[]`/`[1,2,2,3]`/`[2]`; `analysis_quality()` is still arm-independent. The "redesign" is prose only |
| msprime cross-check rerun | **NO** | `s4_msprime.py` (03:39) / `s4_msprime.json` (03:52) unchanged; numbers identical to round 1 (0.716/0.627); the script still scores negatives from the attacker's own reference panel (line 76: `af_lrt_attack(ref["G"], ...)`) |
| F1 dropped, S5 boundary added | **YES** | new `s5_boundary.py`, disjoint negatives, clean grid |

The doc-code mismatch on S3 is an **integrity item**, not a nitpick: the
results document asserts a redesign that does not exist in the artifacts. Fix
the text or the code before anything ships.

---

## 1. Re-verdict C1/C2/C3

**C1 (credibility) — SUPPORTED** (with stated caveats).
The disjoint-negative fix is real and did exactly what the round-1 review
predicted: AUROCs roughly halved (F5 0.92→0.72, F3 0.91→0.71), the tier
ordering survived (F4 0.517 < F3 0.712 ≈ F5 0.720 < F5_sitewise 0.859), the
ancestry-only control stays at chance (0.495 at every tier), and TPR@1% CIs
exclude zero. The honest restatement — "the earlier numbers were inflated by
~0.2 AUROC" — is exactly the accounting a defensible paper needs. Caveats to
carry: **2 seeds** (down from 5; the CI resamples individuals within pooled
scores and does not capture seed-to-seed variability), per-ancestor breakdown
absent from S1b, and the Beacon family remains dead-by-construction (unchanged
strawman — one working attack family carries the claim).

**C2 (spine) — SUPPORTED**, with a mandatory mechanism correction.
+0.108 is real: 3/3 seeds (0.850/0.861/0.871 vs base 0.749/0.749/0.759 — every
+M2 seed exceeds every base seed), and it *grew* after the fix because the
honest base (0.753) no longer masks it. λ_GC/condition-bit at exactly 0 is
correct (constants). But my diagnostic shows the mechanism is **not** "the PCA
spectrum carries extra population information" — it is **in-sample-ness of the
study-fitted spectrum** (details in §3). The claim as worded survives; the
paper's interpretation section must be rewritten around the actual mechanism,
which is both more precise and, frankly, a better story (§4).

**C3 (bet) — NOT SUPPORTED** (unchanged from round 1).
No new S3 evidence exists. Probe counts remain hardcoded inputs, quality
remains arm-independent by construction, `cumulative_auroc` remains a max over
steps, and the sanitized arms' AUROCs still re-state S1's precision tiers. The
results doc's new "coarse-first search, equal quality" wording describes the
same hardcoded behavior with better prose. As submitted, C3 has zero empirical
content — and now additionally carries the doc-code mismatch risk above.

---

## 2. Did the disjoint-negative fix change confidence, and in which direction?

**Up, substantially — for C1 and C2.** Three reasons:
1. The predicted signature appeared precisely: absolute AUROCs fell ~0.2, the
   ordering and the chance-level ancestry control survived. The round-1
   numbers are now *explained* (quantified bias) rather than merely doubted.
2. The fix was implemented correctly (verified in code, and the S1b/S2b
   outputs are consistent with it).
3. Removing the bias *strengthened* C2: the Δ grew from +0.037 to +0.108 once
   the base was de-inflated — the PCA-spectrum channel was being compressed by
   the ceiling effect of the inflated base.

Confidence in C3 is unchanged (no new evidence). One caveat on the fix itself:
2 seeds is thin for a "primary result" table; restore 5 seeds before
submission (minutes of CPU — the round-1 sweep already ran 5).

---

## 3. Is Δ+0.108 credible, or is a residual artifact present?

**It is a genuine, attacker-implementable signal — but its mechanism is
in-sample-ness of the study-fitted spectrum, and the paper must say so.**
Evidence from the diagnostic (3 seeds, exact reproduction of S2b):

- **T1 (ablation):** `pcdist` (projection distance from the study centroid)
  *alone* gives AUROC **0.79–0.81** — stronger than the final AF release
  (0.753). `pcalign` (|PC1| alignment) alone is at chance (0.49–0.51). The
  entire Δ is the PC-distance channel.
- **T2 (independent-spectrum control):** refit the loadings on the attacker's
  reference panel (same size, same population; both target panels
  out-of-sample) → base+M2 drops to 0.748–0.763, **Δ ≈ 0** (−0.001 to +0.014).
  The +0.108 exists *only* when the spectrum was fit on the study members.
- **T3 (attacker-realistic centering):** centering both panels by the public
  study mean (what a real attacker would do) leaves Δ intact (0.102–0.112) —
  the attack is implementable from the release alone.
- **T5 (sanity):** within-ancestry label permutation → 0.49–0.52.
- **T4 (mechanism):** members sit systematically *farther* from the study
  centroid than non-members (mean 0.1448 vs 0.1332; per-individual effect
  ≈ 0.8 sd; panel means differ by ~25 standard errors). Σ top-10 eigenvalues
  of the study sample vs an independent sample are equal (ratio 0.999), so
  this is not sample-eigenvalue inflation — it is **direction tuning**: the
  study's fitted PCs maximize the member sample's variance along them
  (in-sample explained variance = λ̂ exactly), so any independent
  same-population panel projects with systematically smaller spread.

**Consequences for the paper:**
1. Reframe C2: "a study-fitted PCA spectrum is membership-testable *by itself*
   — the projection-distance attack (0.79) outperforms the released AF-LRT
   (0.753), and combined they reach 0.86." That is a stronger governance
   finding than "+0.108 marginal".
2. Report the ablation (pcdist vs pcalign) — pcalign is inert.
3. Report the ref-fit control as the mechanism identification.
4. Flag the magnitude as dimension-dependent (m/n ≈ 2900 SNPs / 900 members
   here): add a cheap m/n sweep (vary n_snps and n_per_group with everything
   else fixed) so the claim is not anchored to one operating point.
5. The direction ("members are farther from the study centroid") should be
   stated and interpreted — it is counterintuitive and reviewers will ask.

---

## 4. Ranking by reviewer value

1. **C2** — now the paper's headline, *if* the mechanism section is written.
   Clean 3/3-seed result, attacker-implementable (T3), mechanism identified
   (T2), and the reframed story ("the diagnostics release is the largest leak,
   stronger than the summary release") is novel, quantified, and directly
   answers the governance question the paper exists for.
2. **C1** — now defensible and honestly accounted (inflation quantified at
   ~0.2 AUROC); the site-size × precision boundary (S5) is the actionable
   artifact. Still a confirmation of 2008-era results, so it is the spine, not
   the contribution.
3. **C3** — zero evidential weight, plus an active doc-code mismatch. Shipping
   it as-is would cost credibility with the same reviewer who checked S1b.

---

## 5. Single most important remaining experiment; strongest surviving rejection

**Most important remaining experiment: S3, actually implemented — or excised.**
Either build the analyst the results doc already describes (a stopping rule
that is *measured*, probe behavior emergent, quality a function of the reveal
policy — an LLM-driven or utility-maximizing analyst, ≥2 families × ≥3 seeds),
or delete C3 from the paper and the "non-degenerate S3" sentence from the
results doc. There is no third option: the current text claims work the
artifacts do not contain. Two cheap must-dos alongside it, both ablations of
what exists rather than new modules: (a) promote the diagnostic's ref-fit
control + pcdist/pcalign ablation into the paper's S2b mechanism section
(`results/raw/_s2b_diagnostic.py`, minutes of CPU); (b) re-run S4 with
`split_population3` — its current "directionally consistent" cross-check
compares a biased instrument (ref-panel negatives) to an unbiased S1b, and the
doc's "rerun" claim is false as written.

**Strongest surviving rejection argument:** "The paper's third pillar (C3, the
governance bet) is a scripted scenario whose probe counts are hard-coded
inputs — and the results document describes a redesign that does not exist in
the code. Its msprime cross-check still uses the self-referential control the
authors themselves flagged. What remains is one clean negative-control result
(C1, confirming 15-year-old leakage results) and one novel quantitative claim
(C2) whose mechanism — in-sample-ness of the fitted spectrum — the authors
apparently have not identified, and whose magnitude is an artifact of their
simulator's dimensionality. All evidence is simulated; the second attack
family is non-functional by construction."

**Cheapest pre-emptions (no new modules):** fix the doc-code mismatch (text
edit), re-run S4 with the disjoint split (minutes), restore 5 seeds in S1b
(minutes), add the mechanism/ablation subsection to S2b (the diagnostic is
written), and add the m/n sweep (minutes). Each removes one clause of the
rejection above. The only clause that requires *new* experimental work is the
S3 analyst — hence it is the single most important remaining experiment.

---

## Minor items carried from round 1 (unfixed, harmless but fix in passing)

- S0 `min_p_causal: NaN` still unfootnoted (why S3's "quality = 5" counts 5 of
  10 planted causal SNPs).
- S2b docstring promises "permutation decomposition and bootstrap CIs"; output
  has neither. Deliver or reword.
- S1b ancestry control compares ref-panel vs non-member-panel distances (both
  non-member panels) rather than members vs non-members; valid as a panel-
  identity control, but say so in the caption.
- `stratified_permutation_auroc` in common.py remains unused by any output
  (a round-1 plan pre-registration still unexecuted).
- Beacon family unchanged since round 1: common-SNP queries, beacon n=900,
  genome-mean f → p₀≈1 → null by construction. Drop it or implement it
  properly; "structurally weak" still reads as a strawman.
