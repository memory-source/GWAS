# Cross-model review verdict — Audit-Governed Multi-Site GWAS (2026-10-02)

Reviewer scope: hostile-but-fair, PoPETs/USENIX/CCS audience. I verified every
headline number against `results/raw/*.json` (all match the bundle and
`EXPERIMENT_RESULTS.md`) and read all five scripts in `code/` (common.py,
s0–s4) plus `refine-logs/EXPERIMENT_PLAN.md` against which the run should be
judged. Verdicts below are grounded in the code, not the prose.

---

## 1. Verdict on C1–C3

**C1 (credibility) — PARTIAL.**
Direction is supported: released AF leaks individual membership with the
ancestry-only control at chance (0.492 everywhere, 5 seeds, per-ancestor
uniform 0.92–0.95). But the headline magnitudes are inflated by design:
`attack_af()` uses the attacker's *own reference panel* as the non-members
(`sn = af_lrt_attack(ref["G"], ...)` where `r_af = allele_freq(ref["G"])`).
Ref-panel individuals pulled toward `ref_af` get systematically *lower* scores,
so the test measures "which of two panels does this individual belong to"
(mean separation δ_member + δ_ref ≈ 2δ) rather than the canonical one-sided
membership test the claim names (separation δ alone, non-members disjoint from
both panels). Additionally the second attack family (Beacon, 0.511) is null *by
construction* (see §5/§6), so one family carries the entire claim. Existence of
AF-based membership leakage is 2008-textbook (Homer; Sankararaman); C1's value
rests entirely on the frontier quantification, which inherits the bias.

**C2 (spine) — PARTIAL.**
The +0.037 ΔAUROC is statistically solid *within this setup*: every
+PCA-spectrum seed (0.950/0.942/0.954) exceeds every base seed
(0.921/0.900/0.914), seed-std 0.005, direction consistent; the CV protocol is
clean (StratifiedKFold, scaler fit on train folds only, out-of-fold
predictions). Reporting λ_GC/condition-bit as exactly 0 individual leakage is
correct — though note a constant feature *cannot* change AUROC, so the "honest
zero" is trivially true by construction, not an empirical null. Two gaps: the
same self-referential negative set applies, and no ablation identifies which
M2 component (pcdist vs pcalign) drives the increment. The docstring promises
"bootstrap CI"; none is in the output.

**C3 (bet) — NOT SUPPORTED** as empirical evidence.
The probe counts (3.0/0.0, 4.0/1.0) are **hardcoded constants** in
`probe_sequence()` — the dependent variable is an input. `analysis_quality()`
is arm-independent (returns causal∩QC-mask count; its own docstring says
"equal quality by construction"), so the "at equal quality" clause is vacuous.
`cumulative_auroc` is a *max over probe steps*, not a cumulative measure. The
noisy-refinement branch never fires in the reported seeds (probe counts are
exactly the base sequences). The sanitized arms' AUROCs merely reproduce S1's
precision tiers (A-sanitized 0.556 ≈ F4's 0.562; B-sanitized 0.911 ≈ F3's
0.911) — i.e., S3's only measured content re-states S1. The "pre-registered
conviction rule" (≥20% reduction) is guaranteed to pass when probe counts are
hardcoded; pre-registering a constant is not a conviction rule. The doc's
"scripted analyst" caveat understates the situation: the effect is *assumed*,
not measured.

---

## 2. Single most important missing experiment

**Re-run S1/S2 with a third, fully disjoint non-member panel** (negatives drawn
from neither the study nor the attacker's reference). The EXPERIMENT_PLAN
pre-registered "attacker reference panel disjoint vs overlapping (report
both)"; only the self-referential condition was run. The simulator already
generates 800 unused individuals per group (`build()` uses n_per_group=1400,
the split consumes 600), so this is a control-group fix at near-zero compute
on the same seeds. Every headline AUROC (0.91–0.99) is currently measured on
the biased instrument; this one change determines whether the frontier's
*absolute levels* are citable. My expectation: the tier *ordering*
(F4 < F3 < F5 < F5_sitewise) and the chance-level ancestry control survive,
absolute AUROCs drop meaningfully — and C1/C2 then stand on defensible numbers.

Runner-up (external validity, next round): one real-data demonstration, e.g.
1000 Genomes super-populations as pseudo-sites attacking a real released AF
vector, to blunt the all-simulated rejection. Fix internal validity first.

---

## 3. Sanitized-diagnostics paradox: real finding or model artifact?

**Model artifact — more precisely, an assumption instantiated as constants.**
Three independent pieces of evidence: (a) probe sequences are hardcoded, so the
100%/75% reduction is an input, not an output; (b) quality is arm-independent
by construction, so "equal quality at fewer probes" is tautological; (c) the
only measured quantities (sanitized-arm AUROCs) are S1's F4/F3 tier attacks
re-run. The idea itself is plausible and unrefuted — pre-providing sanitized
diagnostics could plausibly satisfy an analyst's information need and reduce
probe reveals — but as it stands S3 demonstrates "if the analyst stops when we
assume it stops, final-release exposure equals the precision tier of that
release," which is a corollary of S1, not a behavioral finding. To make it
real: an analyst whose stopping rule is *measured* (utility-maximizing or
LLM-driven, probe behavior emergent) with quality measured as a function of
the reveal policy, not fixed. Note also that even taking S3 at face value,
family B's exposure reduction is ~0.007 AUROC (0.928→0.921) despite 75% fewer
probes, and family A's drop (0.928→0.556) is driven by capping the release at
AF-1dp — i.e., the exposure story is carried by precision capping (an S1
result), not by probe-count reduction per se.

---

## 4. Ranking by reviewer value

1. **C2** — the only genuinely novel, policy-relevant quantitative claim:
   diagnostics reveals carry measurable marginal membership leakage beyond the
   final summary release (+0.037), with study-level scalars honestly at zero.
   This is the question a reveal-governance paper exists to answer. Value is
   contingent on the disjoint-negative fix and one ablation of the M2 driver.
2. **C1** — the cleanest-looking design (ancestry-matched, per-ancestor
   uniform, site-size boundary is a nice artifact: TPR@1% 0.89→0.27 across
   n=50→900), but the phenomenon is established and the control is currently
   self-referential. After the fix it is the solid, unsurprising spine.
3. **C3** — most novel-sounding, currently zero evidential weight. Would rank
   first if the analyst behavior were measured; as submitted it is the weakest.

Flip risk: if the disjoint-negative re-run collapses the +0.037, C1's frontier
becomes the paper's main asset and C2 demotes to a negative/documented result.

---

## 5. Strongest rejection argument + cheapest pre-empting experiment

**Strongest rejection:** "The paper's only original quantitative content — the
release-tier frontier — is measured on a self-referential control: non-members
are drawn from the attacker's own reference panel, converting the one-sided
membership test into an easier two-panel distinguishability test (≈2× the mean
separation). Since AF-based membership leakage is established (Homer'08,
Sankararaman'09), the second attack family is non-functional by construction
(common-SNP queries against a 900-individual beacon with genome-mean AF give
p₀≈1), and S3's probe counts are hardcoded, the paper reduces to one small Δ
(+0.037) measured on the same biased instrument, plus a scripted anecdote."

**Cheapest pre-empting experiment:** the disjoint-negative re-run of S1/S2
(third panel, same seeds, minutes of CPU). If the tier ordering survives and
the ancestry control stays at chance, the rejection's factual core dissolves
and C1/C2 survive with defensible magnitudes. Do not add modules — no MPC, no
real cohort needed for this round; fix the control group on what exists.

---

## 6. Minor hygiene notes (fix in passing)

- "TPR@1% falls monotonically" is false at the tail: 600→900 goes
  0.267→0.274. Say "non-increasing up to n=600" or drop "monotonically".
- The site-size sweep uses 3 seeds; the tier table 5. State both.
- S2 docstring promises bootstrap CI; output has none. Add it or drop the claim.
- Beacon family: implemented with MAF≥0.05-filtered common SNPs, beacon n=900,
  single genome-mean f → p₀≈1 → null by construction. Either implement
  properly (rare-variant beacon, small n, per-SNP LRT) or drop the family;
  "deliberately the weaker one" reads as a strawman to this audience.
- Plan pre-registrations partially unexecuted: stratified-permutation AUROC
  decomposition (`stratified_permutation_auroc` exists in common.py, unused in
  any output), S1 TPR@FPR≤0.1%, and the disjoint/overlapping reference
  comparison. Either run them or say why not.
- S0 `min_p_causal: NaN` — a planted causal SNP fell below the MAF floor.
  Harmless, but it is why S3's "quality = 5" counts only 5 of 10 causal SNPs;
  footnote it.
- S3 uses `n_per_group=800`/split 300+300 while S1 uses 1400/300+300; the
  unsanitized 0.928 vs S1's F5 0.916 are therefore not directly comparable —
  pin the configs in the paper's tables.
