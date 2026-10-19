# Code review verdict — S0–S4 experiment implementation (round: disjoint-negative re-deployment)

Reviewed: `code/common.py`, `code/s0_sanity.py`, `code/s1c_frontier.py`, `code/s2c_mechanism.py`,
`code/s3b_paradox.py`, `code/s4_msprime.py` against the condensed plan. All findings below were
verified by executing the relevant code paths (not just reading).

## Verdict: 1 CRITICAL, 1 MAJOR, 4 MINOR. Checks 1, 2, 3, 5 otherwise PASS.

---

## PASS — the four other review checks

1. **Disjoint-negative split (check 1) — CORRECT.** `split_population3` slices
   `cand[0:300] / [300:600] / [600:900]` per ancestry from `np.where(anc==k)` (ascending), so
   site/ref/nonmem are pairwise disjoint (verified: 0 pairwise index intersections; sizes
   900/900/900 at `n_per_group=1000`). S1c/S2c/S4 all score the **third** non-member panel:
   `s1c` `sn = af_lrt_attack(nonmem["G"], …)`, `s2c` `base(nonmem)` with
   `y = [ones(site), zeros(nonmem)]`, `s4` `sn = af_lrt_attack(nonmem["G"], …)` after its own
   disjoint `idx[:100]/[100:200]/[200:300]` split per msprime population. The attacker's `ref`
   panel contributes only `af_r`/loadings. S0 still scores ref-as-negatives, but the plan only
   requires the disjoint fix for S1c/S2c/S4, and S0 is a meter-sanity script.
2. **AF-LRT sign/orientation (check 2) — CORRECT.** Score = Σⱼ[log P(g|study_af) − log P(g|ref_af)]
   with correct HWE probabilities (p0=(1−f)², p1=2f(1−f), p2=f²); members contributed to
   `study_af` → members score higher. All call sites pass study AF first, ref AF second; the
   rounded tiers (F4 1dp / F3 2dp) don't flip orientation.
3. **AUROC / TPR@FPR / CV (check 3) — CORRECT.** `tpr_at_fpr` maximizes TPR subject to
   FPR ≤ target over unique thresholds. S2c CV is leak-free: `StratifiedKFold(shuffle=True,
   random_state=seed)`; `StandardScaler` and `LogisticRegression` are fit on `X[tr]` only;
   out-of-fold predictions are scored. Bootstrap CI resamples pos/neg within labels.
   (One caveat on tie handling → MINOR-3.)
4. **Seeds & JSON (check 5) — CORRECT.** s0 seed 0; s1c/s2c seeds 0–4 (CV seeded per seed);
   s3b/s4 seeds 0–2; msprime uses `seed+1`/`seed+2`; bootstrap fixed at seed 0. All results
   written to `results/raw/*.json`, and the paper's Table 1 / §mechanism numbers match the JSONs.

Also verified: S3b is genuinely non-degenerate (stopping rule computed from the revealed
spectrum; quality = re-run association at chosen k; no hardcoded probe lists), and the S2c
control logic T1/T2/T3 is implemented as described (T2 ref-fitted loadings collapse Δ to ≈0;
T3 public-mean centering leaves Δ intact — matches the paper's +0.003 / +0.10 claims).

---

## CRITICAL — S0/S1c "ancestry-only control" is a no-op; the paper reports its ≈0.49 as a meaningful negative control

**Bug 1 — the study centroid never enters the computation.** `common.ancestry_only_attack`
computes `mu = ref_scores.mean(0)`, but `pca()` centers the genotype matrix, so *every* panel's
PCA scores have mean ≈ 1e−16 (verified: max |mu| = 9.2e−16). The function therefore always
returns −‖target's own centered PCA scores‖ regardless of its second argument — "distance to the
study PC centroid" (paper, §controlled) is never computed.

**Bug 2 — s1c scores the wrong panels.** In `s1c_frontier.py`:

```python
av_m = ancestry_only_attack(release_objects(ref)["scores"], rel["scores"])     # REF panel scored
av_n = ancestry_only_attack(release_objects(nonmem)["scores"], rel["scores"])  # NONMEM panel scored
anc.append(auroc(np.concatenate([av_m, av_n]), lab))   # lab = ones(site members), zeros(nonmem)
```

`lab` assigns label 1 to *study members*, but position-wise it is applied to **ref-panel**
scores; label 0 (non-members) is applied to nonmem scores. Both scored panels are non-members;
study members are never scored. The control therefore measures sampling noise between two
same-population non-member panels → ≈0.5 **by construction** (reported 0.494 at every tier).
S0 has the same defect inverted: the label-1 block contains ref-panel scores and the label-0
block site-panel scores (reported 0.4993).

**Why it matters.** The paper (main.tex, Table 1 caption) states *"The ancestry-only control is
at chance (≈0.49) at every tier"* and §controlled calls it a *"confound upper bound"*. The honest
control — project both target panels onto the **study-fitted** loadings, center by the public
study mean, score by distance to the study centroid — is exactly the S2c `pcdist` feature, and it
is **not** at chance. Measured at s1c scale (seeds 0–2): AUROC ≈ 0.14 with members as positives
(≈0.86 sign-flipped), because members sit systematically **farther** from the study centroid
(in-sample spread inflation — the paper's own §mechanism finding, pcdist 0.791 / oriented 0.86).
So the honest PC control *exceeds* the AF-summary attack (0.72) and matches F5_sitewise (0.858);
only the meta-attack (0.86) is at/above it. The caption claim is contradicted by the paper's own
S2c analysis, and the research-contract C1 ordering ("control weaker than the real attack") fails
on most tiers once the control is implemented as described.

**Exact fix.**

In `common.py`, make the centroid explicit:

```python
def ancestry_only_attack(target_proj, mu):
    """Score targets by negative distance to an explicit centroid `mu`
    (caller projects raw genotypes onto the study-fitted loadings first)."""
    return -np.linalg.norm(np.atleast_2d(target_proj) - np.asarray(mu)[None, :], axis=1)
```

In `s1c_frontier.py` (inside the per-seed tier loop):

```python
Gs = site["G"].astype(np.float64); Gn = nonmem["G"].astype(np.float64)
mu_study = Gs.mean(0)
_, scores_s, Vt_s = pca(Gs, k=10)
av_m = -np.linalg.norm(scores_s, axis=1)                    # members' released PC coordinates
av_n = -np.linalg.norm((Gn - mu_study) @ Vt_s.T, axis=1)    # non-members projected on study PCs
anc.append(auroc(np.concatenate([av_m, av_n]), lab))
```

Same construction for S0 (members in the label-1 block, ref projected on the study basis in the
label-0 block). Then re-run S0/S1c and update the Table 1 caption and §controlled text: the
expected honest value is ≈0.14/0.86, i.e. the PC channel *dominates* the AF channel at most
tiers — consistent with S2c — and the "confound at chance" framing must be replaced.

---

## MAJOR — `maf()`/QC silently drops every SNP with alt-freq > 0.475

`common.py`:

```python
def maf(G):
    return np.minimum(G.mean(0), 1 - G.mean(0)) / 2.0 * 2.0  # minor allele freq
```

`G.mean(0)` is **2f** (mean genotype), not an allele frequency, so this equals min(2f, 1−2f) and
`qc_mask(maf_floor=0.05)` keeps f ∈ [0.025, 0.475] — dropping every SNP with f > 0.475 (true MAF
up to 0.95 should be kept) and keeping sub-0.05-MAF variants down to f=0.025. The same buggy
formula is inlined twice in `s4_msprime.py` (`np.minimum(G.mean(0), 1 - G.mean(0)) / 2`).
Verified impact: S0 reports `n_snps_tested = 1874/4000` (~47% kept, matching P(f ≤ 0.475) for
f ~ U(0.05, 0.95)).

The mask is applied symmetrically (study AF, ref AF, attack), so AUROC *direction* is not biased
— but every absolute number (AUROC, TPR@1%, λ_GC, min-p, planted-signal recovery) is measured on
a release that omits ~half the common variants, which is not the "standard QC → PCA →
association" pipeline the paper describes, and it removes ~half the planted causal SNPs from
testing (S0's `causal_in_mask` conditions on the survivors).

**Exact fix.** In `common.py`:

```python
def maf(G):
    f = G.mean(0) / 2.0
    return np.minimum(f, 1.0 - f)
```

and in both `s4_msprime.py` sites replace `np.minimum(G.mean(0), 1 - G.mean(0)) / 2` with
`np.minimum(G.mean(0)/2, 1 - G.mean(0)/2)`. Then re-run S0–S4 (numbers will shift; the
disjoint-negative comparison must be re-verified on the corrected mask).

---

## MINOR

1. **s2c docstring inverted (code comment only; paper is correct).**
   `s2c_mechanism.py` says members' projections "sit systematically **closer** to the released
   centroid". Measured: members are **farther** (mean pcdist −0.1448 vs −0.1332; in-sample
   spread inflation), and the paper correctly says the independent panel "projects with
   systematically smaller spread". The CV AUROC is sign-agnostic so all s2c numbers stand; fix
   the docstring to "farther / larger spread".
2. **Dead computation in `af_lrt_attack` (common.py).** The first
   `ll_study = _hwe_geno_loglik(g[:, None, :], study_af[None, :])` allocates an (n, 1, m) array
   that is immediately overwritten by the correct line below it. Delete it (wasted memory/time).
3. **AUROC tie handling.** `auroc` assigns ordinal (not mid) ranks via stable argsort, so
   pos/neg ties count as losses. Harmless for these continuous scores; use midranks
   (`scipy.stats.rankdata`) if scores could ever be discrete.
4. **S3b "100% probe reduction" is by construction, not emergent.** The sanitized arm
   unconditionally returns 0 probes (`probe_until_stable` returns `(0, eig)` when
   `has_sanitized`), so the `both_families_ge_20pct` gate is tautological. Fine as a design
   (pre-provided sanitized block), but don't present the 100% reduction as a measured effect.
   Also `cumulative_attack` uses `mask=None` (no QC mask) unlike every other script, and the
   double `np.round(np.round(af, d), d)` is redundant.

---

## Bottom line

- Checks 1, 2, 3, 5 pass: the disjoint-negative fix is real and correctly wired into S1c/S2c/S4,
  the AF-LRT orientation is right, metrics and the S2c CV are sound, seeds are fixed and results
  are saved as JSON.
- **Do not re-deploy yet.** The ancestry-only control in S0/S1c is a no-op that never scores
  study members and never uses the study centroid; its ≈0.49 feeds a paper caption ("control at
  chance at every tier") that the honest control (≈0.86 oriented) and the paper's own S2c
  mechanism analysis contradict. Fix CRITICAL + MAJOR, re-run, and revise the Table 1 caption /
  §controlled narrative before submission.
