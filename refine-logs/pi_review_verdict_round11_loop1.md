# Round-11 Review — Loop 1 (run_20261006_f68b0cd6)

**Reviewer**: pi (senior-ML-reviewer persona, artifact-level verification; every number below was
checked directly against `results/raw/*.json` and `code/*.py`, not against executor notes)
**Scope**: round-11 limitation-closure revision (s37–s42, enforce.py, genetic_map.py, semi_synthetic.py,
merged s37 shards) + regression check of the paper text built on rounds 8–10.
**Date**: 2026-10-06

---

## 1. What was verified and HOLDS

All six round-11 artifact groups exist, ran, and their headline numbers match the paper to the
quoted precision (verified value-by-value against the JSON keys):

| Paper claim | JSON key | Paper | Artifact | Match |
|---|---|---|---|---|
| s39 library bypass lands on disk | `library_leak_frac` (3 seeds) | 3/3 | 1.0 | ✓ |
| s39 isolated natural-egress blocked / sanctioned push blocked / declared flow / chain | 4 frac keys | 3/3 each | all 1.0 | ✓ |
| s39 transported governed tier | `governed_af_auc` / `governed_pc_auc` | 0.624 / 0.509 | 0.6239 / 0.5094 | ✓ |
| s38 loadings A*, ε_lb(.01) | `pc_loadings.ks_advantage` / `eps_lb_0.01` | 0.61 / ≈24 | 0.6102 / 23.93 | ✓ |
| s38 af_1dp A*, ε_lb | `af_1dp.ks_advantage` / `eps_lb_0.01` | 0.19 / 1 | 0.1871 / 1.057 | ✓ |
| s38 controls | both controls | A*≤0.04, ε≤0.24 | 0.0433/0.2438, 0.0327/0.1915 | ✓ (2dp) |
| s40 bundle vs max single, excess, 5/5 | `bundle`/`max_single`/`excess`/`excess_frac_gt0` | 0.962/0.947/+0.015/5/5 | 0.9619/0.9471/0.0149/1.0 | ✓ |
| s40 LOO per-site-AF cost, noise control | `loo.af_site`, `af_plus_noise` | 0.076 / 0.839 | 0.8863 (cost 0.0756) / 0.8389 | ✓ |
| s41 planted 10 cM on chr1/2/22 | `planted_lsr_auc` | 0.715/0.729/0.702 | 0.7149/0.7294/0.7016 | ✓ |
| s41 real arm ≤ AF | `lsr_auc` vs `af` | 0.456/0.479/0.424 vs 0.866/0.864/0.838 | exact | ✓ |
| s41 50k markers, 3 splits, ~10 cM plant | `n_markers`, `n_splits`, `planted_seg_cM` | 50k/3/10 cM | 50000/3/9.98–10.04 | ✓ |
| s42 ρ sweep + ρ=0 null | `grid.rho` | 0.502/0.641/0.762/0.891/1.000 | exact | ✓ |
| s42 bonus sweep | `grid.bonus` | 0.50→0.58→0.85→1.00 | 0.5016/0.5801/0.8478/1.0 | ✓ |
| s37 gap −0.03…−0.10, AF first at all 6 cells | `delta_pc_minus_af` | −0.03 to −0.10 | −0.0296…−0.0958 | ✓ |
| s37 in-sample gap 0.20→0.08→0.03→0.00 | `pc_in_sample_gap` | same | 0.1989/0.0813/0.0250/0.0002 | ✓ |
| s37 relatedness 5/6 within ~10%, nsp150 +11% | `relatedness_*` | same | 141.2 vs 127.7 = +10.6% | ✓ |
| s37 LD / AF-spectrum / eig10 controls match donor | `controls` | match | within 2–4% | ✓ |
| Genetic map is real | `data/external/genetic_map/*.txt` | IMPUTE2 map | 43k–279k rows, real format | ✓ |
| s37 caches back n_seeds_used=3 | `results/cache/semisyn_*` | 3 seeds/cell | 3 npz per delivered cell | ✓ |
| Stale "remains argued, not built" limitation | paper text | removed | grep = 0 hits; "built and tested" present | ✓ |
| Compile 15 pp, figures present, no active TODO | main.log / figures/ | 15 pp | 15 pages, 8 PDFs, 0 \todo | ✓ |

Regression spot-checks of the older-round numbers quoted in the same sections also matched:
s23 (0.845/0.882→0.622/0.510), s24 (0.941 vs 0.882, 5/5), s25 (ΔE +0.24 CI [0.23,0.26]; null family
−0.24 at p_san≥2; `null_family_paradox: false` at the operating point), s26 (0.866/0.486/0.463/0.653),
s21 (monotone at all 3 n; crossing ≈0.043–0.052 measured; realized 0.05 = 0.63× nominal; real ≈0.10),
s22 (−0.115→+0.096), s10 (0.879±0.009 vs 0.832±0.021, +0.046 CI, +0.103 CI, 25/25),
s12 (+0.163±0.036 at (400,500); −0.012±0.011 at (1600,8000)), s16 (0.768/0.675/0.805; 0.727/0.514/0.512),
s17 (0.909/0.896; 0.807/0.745; +0.056), s18, s19 (0.990/0.487), s30 (0.66–0.67 flat; 2.2k→48k),
s31 (0.708/0.611; 0.780/0.611), s36 (0.671/0.486), s15/s20 (fst-like 0.0994/0.0503; 0.1418/0.1139;
top-10 var 0.258/0.142). s29 slow cell: −0.0136, P=0.21, 200 seeds ✓.

**The five round-11 limitation-closure claims are real and supported by their artifacts.**
The enforcement bypass is genuinely demonstrated then closed; the bundle has the right
falsification control; the IBD bound uses a real genetic map with a planted positive control on
every chromosome; the missingness curve has a ρ=0 null exactly at chance; the semi-synthetic
cohort assigns labels post-generation (no circularity) and carries the four mandatory controls.

---

## 2. Score: 7.5 / 10 (current draft); 8.5 attainable after the fix list below

Deductions are text-level: one claim contradicted by the paper's own artifact, one absolute
number whose estimator is unsound as described, and one internal contradiction. The
experimental substance behind every section verified clean.

## 3. Remaining weaknesses (ranked) and MINIMUM fixes

**W1 (moderate — claim contradicts its own artifact).** §governance: "For the leaking channels
$A^\ast$ exceeds the binormal value $2\text{AUC}-1$." The JSON says the opposite:
`gap_vs_2aucm1` is **negative** on all five leaking channels (−0.056…−0.153) and positive only
on the two controls (+0.029, +0.043). `EXPERIMENT_RESULTS.md` B1 states the correct direction
(">0 on the two controls, <0 on the leaking channels"); the paper sentence inherited the plan's
incorrect success criterion (`EXPERIMENT_PLAN.md` B1: "A*−(2AUC−1)>0") and was never checked
against the delivered JSON.
**Minimum fix**: swap the sentence — "A* differs from the binormal value 2AUC−1 (below it for
the leaking channels, above it for the at-chance controls), so the scores are not Gaussian and
A* is not AUROC in disguise."

**W2 (moderate — unsound absolute number).** The ε_lb presentation. (a) The paper's formula
$\varepsilon_{lb}=\max_\alpha\log((\text{TPR}(\alpha)-\delta)/\alpha)$ is the single
TPR-branch only, but the code takes the max of two branches
(`eps_lower_bound`: b1 = log((TPR−δ)/α) **and** b2 = log((1−δ−α)/(1−TPR)), a swapped-hypothesis
branch), and the quoted ≈24 values are **b2-dominated**: they come from the empirical ROC point
where TPR=1 (zero false negatives among 900 positives), where 1−TPR is clipped at 1e-12
(log(0.09/1e-12) ≈ 25). The paper's own formula would give ≈3–4, not ≈24. (b) Calling this "a
finite-sample lower bound" is anti-conservative: with n_pos=900, the true FNR at that operating
point could be ≈3/900 at 95% confidence (rule of three), giving ε_lb ≈ 3, not 24. (c) The stated
n_neg≈300 is wrong — `cell_frontier` scores 900 negatives (300/ancestry × 3).
**Minimum fix** (text-level): state the two-branch formula, fix n_neg to 900, and add one
sentence: the δ=0.01 bound is set by the zero-observed-FN point and should be read as
ordering-only (leaking ≈ tens, 1dp ≈ 1, controls ≈ 0.2), not as a calibrated ε estimate.
Optional (stronger): recompute with a confidence-corrected FNR floor (rule of three).

**W3 (minor-moderate — internal contradiction).** §replication still says "A full
phased-haplotype/IBD attack (Bu et al.) remains future work." while §scale reports exactly that
attack as run (s36 chr22 marker-length; s41 three chromosomes in cM, with numbers). One of the
two must go.
**Minimum fix**: delete the §replication sentence or rewrite it to point forward ("§scale
reports this attack").

**W4 (minor — wrong cross-reference).** Limitations: "The biobank-scale cohort arm (a
semi-synthetic cohort built from real 1000G haplotypes, Section~\ref{sec:governance})" — that
arm lives in §sec:scale, not §governance.
**Minimum fix**: `\ref{sec:governance}` → `\ref{sec:scale}`.

**W5 (minor — oversold sweep).** §scale: "the base missingness level shifts the absolute AUROC
only weakly" — the `base_missing` sweep is entirely at ceiling (0.9992–1.000, as is
`sigma_batch`), so the sweep shows **no** variation and cannot characterize a weak shift.
**Minimum fix**: reword to "at the delivered operating point the channel saturates regardless
of base missingness (0.999–1.000)", or re-run one point below ceiling (e.g. ρ=0.5, bonus 0.01).

**W6 (minor — self-contradictory wording).** §scale: "Extending this … at full marker density,
… (50k markers per chromosome strided uniformly, 3 splits)" — the same sentence contradicts
itself: the parse is capped at 200k markers and then uniformly strided to 50k, which is not
full density.
**Minimum fix**: "at uniform genome-wide density (50k markers per chromosome, uniformly strided)".

**W7 (minor — audit completion claim vs file state).** `EXPERIMENT_AUDIT.md` action item 5 is
marked [x] ("refresh ROUND11_DESIGN.md's stale grid text") and `EXPERIMENT_RESULTS.md` says the
grid text was "corrected to the delivered grid" — but `ROUND11_DESIGN.md` line 29 still
advertises "nsp ∈ {150, 500, 2000, 8000, 33000} × density {thin, dense}", the grid that never
ran (delivered: 6 cells, 150–12000, no thin/dense split). Same doc: B5's "donor coefficient of
variation 0.38 vs cohort 0.32" — the cohort CV is 186.03/141.25 = **1.32**, not 0.32 (donor
0.38 is correct).
**Minimum fix**: actually apply the design-doc refresh; fix the CV figure (or drop the clause).

**W8 (nits).** `s38_privacy_metric.py` imports `anc_only_attack` but never calls it (the audit
removed s40's dead imports but missed this one). `governance.py` docstring still says "The
sole-egress / SPU-binding claim remains argued" — now false (enforce.py builds and tests it).

## 4. Structural caveats (recorded, not new findings)

- Every membership-level arm is simulated or semi-synthetic; fully-real data is 1000G with
  synthetic splits and a constructed missingness channel. The paper discloses all of this
  (semi-synthetic labeled, ρ-parameter honest, IBS arm disclosed as unphased) — the remaining
  exposure is inherent to the measurement design, not a hidden flaw.
- s37/s41/s42 use 3 seeds/split per cell with tight stds on the headline deltas; acceptable for
  the claims made, but a reviewer may ask for more seeds on the s37 relatedness control.
- The s42 "batch-composition" framing is honest (ρ = fraction of member-heavy populations) but
  the channel remains a per-individual call-rate bonus by construction; the ρ=0 null at 0.502 is
  the right control and it holds.

## 5. Verdict

**Almost — not ready as-is; ready after the fix list.** No new experiments are required:
W1–W8 are one-sentence text edits plus one doc refresh (~1–2 hours). After they land, the
paper's round-11 story is fully artifact-backed: the previously argued sole-egress mechanism is
now built and tested (and the stale limitation is already removed), the bundle composition
result has a real falsification control, the haplotype channel is bounded on real physical
length across three chromosomes, and the biobank-scale cohort claim is supported by its
controls. Recommend: apply W1–W8, recompile, then treat round 2 of this loop as a
regression check.
