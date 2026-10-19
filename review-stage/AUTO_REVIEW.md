# AUTO_REVIEW — run_20261004_d429b8c1

Backend: pi-review (local glm-flash), difficulty: medium, MAX_ROUNDS=4.
State: `review-stage/REVIEW_STATE.json`. Traces: `.aris/traces/auto-review-loop/2026-10-04_run01/`.
Reviewer memory: `review-stage/REVIEWER_MEMORY.md`.
Artifacts under review: `paper/main.tex` / `main.pdf` (6 pp),
`code/*.py`, `results/raw/*.json`.

## Round 1 (2026-10-04) — Reviewer backend: pi-review

### Assessment (Summary)
- Score: 7/10
- Verdict: almost
- Key criticisms: (1) release-unit ambiguity — the projection attack needs the
  PC loading matrix, but `release_objects` emitted only eigenvalues and T3
  recomputed loadings from raw genotypes; (2) operating point unstated;
  (3) uncertainty gaps (s7 single split, 2-seed sweeps); (4) text/numeric slips
  ("+0.003±0.006" vs ±0.005; "0.85" vs 0.845; "sign-reversed" wrong; top-k
  criterion undisclosed); (5) F/M taxonomy undefined; (6) missingness row
  near-tautological; (7) "largest at small sites" should be ratios; (8) broken
  `\ref{sec:discussion}`; (9) stale docs.

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response</summary>

# Senior Review — GWAS Leakage Measurement Paper (Loop 2, Round 1)

I read `main.tex`, all 13 scripts + `common.py`, and all 13 raw JSONs in full, recompiled the paper fresh, and script-matched 66 paper values against the JSONs. Verdict file: `refine-logs/pi_review_verdict_round6.md`.

## What holds up (verified, not trusted)
- 64/66 checked paper numbers match `results/raw/*.json` exactly (Tables 1/2/3, per-ancestor, CIs, msprime, real-data, s8, s3b, s1d). Two apparent misses were rounding tolerance.
- Splits are genuinely 3-way disjoint; the s7 1000G split shuffles within super-population before splitting; msprime slices are exchangeable.
- MAF filter is the fixed version; 3814/4000 tested matches "≈3800 of 4000".
- Orientations correct and honestly disclosed; PC naive(0.12)/adversarial(0.88) backed by s0b and s2c.
- Controls real: ref-basis structure control 0.502/0.510, permutation 0.486, T2 refit collapse (+0.003), T3 centering intact (+0.10), eig ratio 1.000, noise 0.496; CV is out-of-fold.
- s3b genuinely implemented; real VCF present; GT parsing fixed. Compiles 6pp/0 errors.

## 1. Score: 7/10
## 2. Critical weaknesses (ranked)
1. "Attacker-implementable from the release alone" rests on an unstated release unit (needs loadings; release_objects emits eig+scores only; T3 recomputes from raw genotypes).
2. Operating point stated nowhere (900/900/900, 4000 SNPs).
3. Uncertainty gaps: s7 single split; boundary/dim 2 seeds unstated; Table 2 no std.
4. Text/numeric slips: "+0.003±0.006" vs ±0.005; "0.85" vs 0.8446; "sign-reversed" wrong (same sign, attenuated); real-data top-200 criterion differs from simulated; Table 2 caption pointer.
5. "Largest intermediate-release channel" comparison class undefined; F/M taxonomy never explained.
6. Missingness row near-tautological (iid uniform noise).
7. "Generally largest at small-to-moderate sites" — pattern cleaner in ratios.
8. Broken `\ref{sec:discussion}` → "Section ??" in PDF.
9. Stale docs (EXPERIMENT_RESULTS.md; tracker "compiled clean" with 3 overfull + 1 undefined ref).
10. auroc() no average-rank tie handling — latent only.

## 3. Minimum fixes
- W1: name the released object (PC loadings) + conditional scoping; optionally add loadings to release_objects and have T3 consume them.
- W2: state the operating point.
- W3: state seed counts; optionally rerun s7 with several split seeds.
- W4/5/7/8: text edits (±0.005; 0.84(5); "attenuated, same direction"; top-k disclosure; define F/M; "ratios"; add label).
- W6: acknowledge structureless missingness.
- W9: sync stale docs.

## 4. READY? Almost

</details>

### Actions Taken
- **W1 (substantive)**: `release_objects` now returns `loadings`; `s2c_mechanism`
  consumes `rel["loadings"]` (the released object) for all study-fitted attacks
  instead of recomputing from raw genotypes; paper scopes the implementability
  claim as conditional on loadings being shared (as federated PCA does).
  Re-ran s2c: numbers unchanged (base 0.839, +PC 0.941, Δ+0.102, T2 +0.003,
  T3 +0.102, perm 0.486, eig 1.000).
- **W2**: added an "Operating point" paragraph (300/ancestry → 900/900/900,
  4000 SNPs ≈3814 tested, ratio ≈13; seed counts per experiment).
- **W3 (substantive)**: re-ran S7 with 5 random three-way splits; real-data
  AUROC now 0.684/0.681/0.564/0.640/0.808 with per-split std; structure control
  0.504. Paper updated (was single split 0.699/0.696/0.577/0.646/0.793).
- **W4**: fixed "+0.003±0.006"→"+0.003", "0.85"→"0.84 at both 2dp and full",
  "sign-reversed"→"same-signed but strongly attenuated" with both orientations.
- **W5**: defined the F/M taxonomy once in §6; noted per-site AF is itself
  intermediate but reported separately.
- **W6**: acknowledged structureless simulated missingness in the Discussion.
- **W7**: reworded dimension-dependence to "ratios" framing.
- **W8**: added `\label{sec:discussion}`; recompiled twice — 0 undefined refs.
- **W9**: synced EXPERIMENT_RESULTS.md (real-data PC leak) and tracker.
- **W10**: noted; auroc tie-handling affects no reported number (continuous scores).

### Results
- Paper: 6 pages, 0 errors, 0 undefined references (recompiled twice).
- S7 rerun (5 splits) and S2c rerun (released loadings) logged in `results/raw/`.

### Status
- Continuing to round 2 (score 7, verdict "almost" → not a stop).
- Difficulty: medium

## Round 2 (2026-10-04) — Reviewer backend: pi-review

### Assessment (Summary)
- Score: 7.5/10 (up from 7.0)
- Verdict: almost (one clause from ready)
- Verified fixed: loadings now the released object and s2c bit-identical; S7 at
  5 splits (0.684/0.681/0.564/0.640/0.808, ctrl 0.504); all text/numeric slips
  gone; compiles 0 errors / 0 undefined refs.
- Remaining: abstract's unconditional "release alone" clause; tracker S7 row;
  condition-bit M-number nit; 3 overfull boxes.

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response</summary>

# Re-Review (Loop 2, Round 2 — confirmation after fixes)

I re-read the three changed code files in full, both re-generated JSONs verbatim, the revised `main.tex`, and both tracker docs; recompiled fresh (pdflatex ×2); and number-matched the new values. Verdict file: `refine-logs/pi_review_verdict_round7.md`.

Fix verification: (a) s2c now consumes rel["loadings"], re-run bit-identical (pcdist 0.8813, ΔT2 +0.0026±0.0047, ΔT3 +0.1020); (b) S7 at 5 splits 0.684/0.681/0.564/0.640/0.808 ±0.009–0.014, control 0.504, PC 0.347/≈0.65; (c) all round-1 slips gone, compile 0 errors/0 undefined refs, 6pp, operating point accurate, F/M taxonomy defined, missingness caveat present.

Score: 7.5/10. Verdict: READY: Almost — one clause from Yes.

Remaining & minimum fixes: (1) abstract add loadings conditional; (2) tracker S7 row one line; (3) nit: condition bit absent from Mn enumeration; (4) cosmetic overfull boxes.

Memory Update: previous suspicions addressed (loadings gap closed, s7 multi-split, slips fixed); no new material suspicions.

</details>

### Actions Taken
- Added the loadings conditional to the abstract.
- Added M4 (condition bit) to the §6 enumeration.
- Synced the tracker S7 row to the new multi-split numbers.
- Recompiled twice: 0 errors, 0 undefined references.

### Results
- Paper: 6 pages, clean compile.

### Status
- **Stopping** — policy-approved positive: score 7.5 ≥ 6 AND verdict "almost".
- Difficulty: medium

## Method Description

The paper is a controlled measurement study of membership leakage from released
statistics in multi-site GWAS. Pipeline: (1) **Panels.** Three ancestries of 300
individuals each are simulated under a Balding--Nichols model ($F_{ST}=0.08$)
and split three ways *disjointly* into study members, an attacker reference
panel, and independent non-member targets; an msprime 3-population demography
and real 1000 Genomes Phase-3 chr22 pseudo-sites (EUR/EAS/AFR) serve as
cross-checks. QC keeps MAF $\in[0.05,0.95]$. (2) **Release objects.** Each run
emits the workflow's released quantities: per-site allele frequencies, the
study-fitted PC basis (eigenvalues, loadings, scores), $\lambda_{GC}$, and a
missingness table, at several precision/visibility tiers (F2 top-$k$; F4 AF-1dp;
F3 AF-2dp; F5 AF-full; F5 per-site). (3) **Attacks.** An allele-frequency
likelihood-ratio membership score plus two functional-form variants
(correlation, learned squared-deviation), and a study-fitted PC
projection-distance attack; an ancestry/structure-only control projects both
target panels onto an *independent* reference basis (chance by construction).
(4) **Measurement.** AUROC and TPR@FPR$\le$1% with bootstrap CIs, per-ancestor
breakdowns, and a cross-validated logistic meta-attack for marginal leakage of
intermediate reveals. Data flow: simulator/real-VCF $\to$ `release_objects` $\to$
attack scores $\to$ metrics JSON in `results/raw/` $\to$ tables/figures in the
paper. All steps are CPU-only, seeded, and traceable 1:1 from paper numbers to
`results/raw/*.json`.

---

## Final Summary (run_20261004_d429b8c1)

Score progression: 7.0 (round 1) → 7.5 (round 2). **Stopped at round 2 on a
policy-approved positive verdict** (score 7.5 ≥ 6 AND verdict "almost").

Result-to-claim verdict (`CLAIMS_FROM_RESULTS.md`): **partial**, confidence high.
The PC-diagnostic dominance claim is supported in simulation (paired 5/5 seeds,
PC 0.881 vs AF 0.845; +0.10 increment; mechanism = in-sample-ness via five
convergent controls); the per-site AF channel (0.964) and the self-reference
bias (0.132) are supported; **PC dominance is not reproduced on real 1000G
chr22** (only the AF-tier ordering is; real PC-distance 0.347 naive / ≈0.65
flipped), and the paper already restricts the PC channel to simulated panels.
The remaining gap to a full "yes" is claim scoping (operating point + attack-side
sign calibration + real-data scope), not new simulation work.

Artifacts: `paper/main.tex` / `main.pdf` (6 pp, 0 errors, 0 undefined refs);
`results/raw/*.json`; `review-stage/REVIEWER_MEMORY.md`;
`.aris/traces/auto-review-loop/2026-10-04_run02/`.

---

# AUTO_REVIEW — run_20261004_bcd63f1a (round 8 re-entry: experiment-gap closure)

Backend: pi-review (local glm-flash, xhigh), difficulty: medium, MAX_ROUNDS=3.
State: `review-stage/REVIEW_STATE.json`. Traces:
`.aris/traces/auto-review-loop/2026-10-04_run03/`. Reviewer memory:
`review-stage/REVIEWER_MEMORY.md`.

Entry context: the previous run stopped at 7.5/"almost" and the result-to-claim
gate returned "partial", naming five evidence gaps. Three sub-agents planned the
next experiments; this run implemented them (s9-s14 plus the s15 control
artifact), updated the paper, then re-entered the review loop at round 1.

## Round 1 (2026-10-04) — Score 7/10, "Almost"

### Assessment (Summary)
- The measured core verifies (≈65 cited numbers match raw JSON; S9c re-run
  bit-identical). Design integrity confirmed (disjoint splits, train-only PC
  calibration, out-of-fold CV, attacks consume released loadings).
- Weaknesses: W1 the real-data "less structured" explanation is refuted by the
  repo's own data; W2 a truncated sentence; W3 the batch-missingness result
  (0.761) was unreported and misread; W4 the real-data paired test was free; W5
  seed-inventory / figure-label items.

### Actions Taken
- W1: replaced both "less structured" sentences with an honest "no simple
  explanation" plus a measured structure comparison (Fst-like 0.10 real vs 0.05
  sim; top-10 variance fraction 0.26 vs 0.14), backed by `code/s15_controls.py`.
- W2: completed the truncated sentence.
- W3: added a Discussion sentence interpreting the batch variant as
  coverage-degraded AF (0.76, below 0.85 full-coverage and 0.80 uniform-70
  controls, `s15_controls.py`), not an additional channel.
- W4: added the real-data paired test (+0.043, CI [+0.033,+0.053], 20/20
  positive, p=1.3e-7).
- W5: seed inventory updated; figure-5(a) labels fixed; "mid-marker" →
  "small-panel".

## Round 2 (2026-10-04) — Score 7.5/10, "Almost"

W1-W4 verified fixed. Remaining: a residual "closer" sign string (should be
members-farther), a "large-panel" regression at the crossover peak, the batch
gloss 0.81 (unbacked) vs the artifact's 0.845, and the unnamed Fst-like
estimator. All fixed in round 3; the estimator is now named and backed by
`s15_controls.py`; the batch sentence cites the three artifact-backed numbers.

## Round 3 (2026-10-04) — Score 8/10, VERDICT: ready

All prior items verified fixed; `s15` re-run from seeds bit-identical. No
blocking weaknesses. One cosmetic item ("mechanism paragraph below"→"above"),
applied. **Loop stopped on a policy-approved positive verdict** (score 8 ≥ 6 AND
verdict "ready").

## Final Summary (run_20261004_bcd63f1a)

Score progression: 7.0 → 7.5 → 8.0; verdict "almost" → "almost" → "ready".

Result-to-claim (round 2, `CLAIMS_FROM_RESULTS.md`): **partial**, high
confidence — a *scoped* partial. All five round-1 gaps closed; the residual is
that the unscoped "PC is THE dominant channel" claim is over-broad: on real
1000G chr22 the calibrated PC channel (0.652) sits below AF (0.688) while adding
a small consistent increment (+0.043, p=1.3e-7); in simulation the PC channel
dominates (paired +0.046, CI [0.037,0.055], 25/25 seeds). The paper's scoped
claim matches what the results support.

Artifacts: `paper/main.tex` / `main.pdf` (8 pp, 0 errors, 0 undefined refs);
`results/raw/*.json` (incl. new s9-s15); `review-stage/REVIEWER_MEMORY.md`;
`.aris/traces/auto-review-loop/2026-10-04_run03/`.

## Method Description

The method is a controlled membership-inference measurement harness for the
*intermediate* objects a federated multi-site GWAS releases (per-site allele
frequencies, the study-fitted PCA spectrum, missingness/QC tables) rather than
only the final association statistic. It builds matched panels — study members,
an attacker reference panel, and an independent non-member set, all drawn
disjointly from the same structured population (parametric Balding-Nichols, or
real 1000G chr22 splits) — computes the released objects with the standard GWAS
pipeline (QC, PCA, marginal association), and scores membership with an
allele-frequency likelihood-ratio attack, a release-fitted PC projection-distance
attack, and a structure-only control that projects both target panels onto an
independent reference basis. Metrics are AUROC and TPR@FPR≤1% with bootstrap
CIs and per-ancestor breakdowns; a cross-validated logistic meta-attack measures
the marginal leakage of each intermediate reveal. Data flow: simulator/real-VCF
→ `release_objects` → attack scores → metrics JSON in `results/raw/` → tables and
figures in the paper. All steps are CPU-only, seeded, and traceable 1:1 from
paper numbers to `results/raw/*.json`.

---

# AUTO_REVIEW — run_20261005_3e9677f2 (limitations closure: bigger data + compute)

Backend: pi-review (glm-flash, xhigh), difficulty medium, MAX_ROUNDS=3.
State: `review-stage/REVIEW_STATE.json`. Traces:
`.aris/traces/auto-review-loop/2026-10-05_run01/`.

Entry context: the paper was at 8/"ready" but its Limitations disclosed (i)
chr22-only real data, (ii) PC channel measured only on the parametric simulator,
(iii) attacks all forms of one marginal channel, (iv) the M3 statistic
structureless, (v) no MPC prototype. This run added data and compute to retire
(i)-(iv): two independent real chromosomes (chr1, chr2), an msprime
recombination-realistic cross-check for the PC channel, a distinct haplotype/LD
attack family, a genuine M3 missingness channel, and a larger site-size sweep.

## Round 1 — Score 8/10, "Almost"
Verified all s16-s20 numbers. Weaknesses: the abstract's ordering claim was
contradicted by its own new numbers (hap-LD 0.99, M3 1.00 > PC 0.88); the s18
"swapped" null was structure-contaminated (msprime populations-as-sites); minor
"+0.056 at both" / docstring / seed-inventory / structure-reconciliation /
operating-point items.

### Actions
- Reran s18's mechanism with an ancestry-matched WITHIN-population split
  (in 0.729 ± 0.015, out-of-sample 0.510 ± 0.010, swapped 0.519 ± 0.012) — a
  genuine de-confounding.
- Re-scoped the ordering claims: "largest among the canonical GWAS diagnostic
  releases" + explicit statement that the haplotype/LD and
  membership-correlated-missingness families can be larger still.
- Fixed +0.055–+0.056, s17 docstring, §4 seed inventory, structure
  reconciliation, hap-LD operating point + release-realism clause.

## Round 2 — Score 8.5/10, "Almost"
Prior fixes verified; two residual unscoped "largest single channel" restatements
(Discussion, §6 intro) — same error class as round 1, caught only by grepping
every instance. Both scoped to "largest marginal-frequency single channel".

## Round 3 — Score 9/10, VERDICT: ready
Full ordering-claim sweep clean; all prior fixes hold; compile 10 pp, 0 errors.
One optional polish (line 412 → "largest canonical intermediate channel"),
applied. **Stopped on a policy-approved positive verdict.**

## Final Summary (run_20261005_3e9677f2)

Score progression: 8.0 → 8.5 → 9.0; verdict "almost" → "almost" → "ready".

New experiments (round 9): `real_multi.py` (chr1/chr2 loader), `s16_multichrom.py`
(chr1/chr2 replication), `s17_scale.py` (site-size sweep), `s18_msprime_pc.py`
(msprime PC cross-check), `s19_haplotype.py` (haplotype/LD channel),
`s20_extra.py` (M3 missingness channel + structure comparison). Results in
`results/raw/s16..s20*.json`.

Key outcomes:
- Real ordering replicates on two independent chromosomes: PC < AF on both
  (chr1 0.675 vs 0.768; chr2 0.689 vs 0.752) but adds a consistent increment;
  the in-sample mechanism and loadings-granularity ladder replicate on all three.
- msprime (recombination-realistic) confirms PC > AF (0.691 vs 0.627) with the
  same granularity ladder and an ancestry-matched in-sample mechanism.
- A distinct haplotype/LD channel reaches 0.990 AUROC but is in-sample-only
  (transfers to a fresh panel at chance, 0.487) — in-sample-ness is a general
  property of fitted released objects.
- M3 missingness is a genuine graded channel when the study batch composition is
  membership-correlated (call-rate excess 0% → 0.49, 1% → 0.93, 2% → 1.00);
  membership-independent missingness is only coverage-degraded AF (0.76).
- New honest boundary: at realistic site sizes (n ≥ 1200, M = 10k) AF saturates
  and PC no longer dominates (AF 0.807 vs PC 0.745 at n=3000); the +PC
  combination increment is stable (+0.055–0.056).

Result-to-claim remains **partial** (scoped): the PC-dominance ordering is a
small-site simulated-regime property; the paper states this boundary explicitly.

Artifacts: `paper/main.tex` / `main.pdf` (10 pp, 0 errors, 0 undefined refs,
Fig. 6 added); `results/raw/*.json`; `review-stage/REVIEWER_MEMORY.md`;
`.aris/traces/auto-review-loop/2026-10-05_run01/`.

## Method Description (updated)

Same controlled membership-inference harness as before, extended in three ways.
(i) Real-data coverage now spans three 1000G chromosomes (chr22, chr1, chr2),
each split three ways disjointly per super-population, with matched operating
points. (ii) The PC channel is additionally characterized on an msprime
three-population demography (real recombination, 5 Mb) and under an
ancestry-matched within-population held-out-site split. (iii) Two further attack
families are added: a haplotype/LD score using the joint distribution of adjacent
markers (fit on the study vs an independent panel), and a genomic-missingness
score using the released per-SNP call-rate table. Data flow unchanged:
simulator/real-VCF → `release_objects` → attack scores → metrics JSON in
`results/raw/` → tables/figures in the paper. All steps CPU-only, seeded, and
traceable 1:1 from paper numbers to `results/raw/*.json`.

---

# AUTO_REVIEW — run_20261005_d6cf307b (governance-layer / model-dependence round)

Backend: pi-review (glm-flash, xhigh), difficulty medium, MAX_ROUNDS=4.
State: `review-stage/REVIEW_STATE.json`. Traces:
`.aris/traces/auto-review-loop/2026-10-05_run02/`.

Entry context: the paper was at 9/"ready" with a Limitations list still naming
(i) model-specific magnitudes, (ii) the unexplained real-vs-sim discrepancy,
(iii) the sanitized-diagnostics paradox as mere conjecture, (iv) no
phased/IBD attack, (v) "the governance mechanism is argued, not built". This run
closed all five with new data and compute: an Fst crossover sweep, a Dirichlet
admixture sweep, a release-object governance prototype (manifest + egress monitor
+ hash-chained ledger) with a violation matrix and a cross-run composition
measurement, the paradox as a pre-registered result with a falsification control,
and an IBS channel.

## Round 1 — Score 8/10, "Almost"
Verified all s21-s26 numbers. Weaknesses: the Dirichlet concentration reading was
inverted; the crossover conflated nominal Fst with measured structure; the "IBD"
arm actually ran on dosages; the anti-composition flag was producer-asserted;
s25 is a corollary of the exposure curve; minor seed/doc items.

### Actions
Reinterpreted the admixture sweep; restated the crossover in measured terms;
reworded the IBS arm; added the ledger-honesty sentence; noted the s25 corollary
and the null-side negative Δ_E; extended the seed inventory.

## Round 2 — Score 8.5/10, "Almost"
5/6 fixed. The crossover conflation survived at summary positions (abstract,
contribution 5, §7 sentence, Fig 7a caption/shade) — the recurring
summary-restatement failure mode.

### Actions
Struck the §7 sentence; fixed the abstract, contribution 5, and the caption;
moved the Fig 7a shade and verified the regenerated PDF.

## Round 3 — Score 9/10, "Almost"
One item: the Fig 7a band was drawn at nominal 0.06-0.10 while the caption
described measured units (the real level, measured 0.10, is nominal ~0.16).

### Actions
Redrew the band at nominal 0.09-0.16 with an explicit "measured 0.06 -> 0.10"
in-figure label; merged redundant §7 sentences.

## Round 4 — Score 9/10, VERDICT: ready
All items verified; two cosmetic one-word polishes applied. **Stopped on a
policy-approved positive verdict.**

## Final Summary (run_20261005_d6cf307b)

Score progression: 8.0 -> 8.5 -> 9.0 -> 9.0; verdict almost -> almost -> almost
-> ready.

New experiments (round 11): `governance.py` (manifest + egress monitor +
hash-chained ledger), `real_phased.py` (keeps VCF phase), `s21_sensitivity.py`
(Fst crossover), `s22_admixture.py` (Dirichlet), `s23_governed.py` (governed vs
ungoverned + injection suite), `s24_crossrun.py` (composition audit),
`s25_paradox.py` (paradox as result), `s26_ibd.py` (IBS channel). Results in
`results/raw/s21..s26*.json`; Fig. 7 added.

Key outcomes:
- The PC-vs-AF ordering is monotone in measured structure and crosses zero at a
  measured statistic ~0.04-0.06, below the real-1000G level (~0.10) — the
  real-data attenuation is the same phenomenon at a different structure value,
  not an unexplained discrepancy.
- A release-object governance prototype blocks/coarsens undeclared or
  over-precision releases (5 blocked + 1 coarsened; chain verifies) and drives
  the leakage attacks to the declared-safe tier; individually declared-safe
  releases compose (cumulative 0.941 > max single 0.882, 5/5 seeds), motivating
  an anti-composition clause.
- The sanitized-diagnostics paradox is now a pre-registered result
  (Δ_E = +0.24 at the precision gap, 0 without it, a null control shows no
  paradox).
- An IBS channel detects planted IBD (0.866) but is at chance unplanted (0.486)
  and on real chr22 (0.463 ≤ AF 0.653), bounding the haplotype-channel claim.

Result-to-claim remains **partial** (scoped): the PC-dominance ordering is a
small-site simulated-regime property; the paper states this boundary explicitly
and now quantifies it as a function of structure.

Artifacts: `paper/main.tex` / `main.pdf` (12 pp, 0 errors, 0 undefined refs,
Fig. 7 added); `results/raw/*.json`; `review-stage/REVIEWER_MEMORY.md`;
`.aris/traces/auto-review-loop/2026-10-05_run02/`.

## Method Description (updated)

The controlled membership-inference harness now additionally (i) sweeps the
simulated structure level (Fst) and a Dirichlet admixture concentration to show
the PC-vs-AF ordering is governed by measured structure, (ii) wraps the
`release_objects()` output in a declarative reveal manifest + egress monitor +
hash-chained ledger and stress-tests it against the standard attack injections,
(iii) measures cross-run release composition, (iv) re-runs the leak attacks on
the governed release, and (v) adds an identity-by-state window-sharing channel
with a planted-IBD positive control and a real-1000G arm. Data flow is unchanged
(simulator/real-VCF -> `release_objects` -> attack scores -> metrics JSON in
`results/raw/` -> tables/figures); all steps CPU-only, seeded, and traceable 1:1
from paper numbers to `results/raw/*.json`.

---

# AUTO_REVIEW — run_20261006_350bbd97 (GPU / scale round)

Backend: pi-review (local glm-flash), difficulty: medium, MAX_ROUNDS=4.
State: `review-stage/REVIEW_STATE.json`. Reviewer memory: `review-stage/REVIEWER_MEMORY.md`.
Artifacts under review: `paper/main.tex` / `main.pdf` (13 pp; new Section
`sec:scale` + Fig. 8), `code/common_gpu.py`, `code/s27..s36*.py`,
`results/raw/s2[7-9]*.json`, `results/raw/s3[0-6]*.json`,
`refine-logs/GPU_ROUND_RESULTS.md`, `refine-logs/EXPERIMENT_PLAN_GPU.md`.

## Round 1 (2026-10-06) — Reviewer backend: pi-review

### Assessment (Summary)
- Score: 7.5/10
- Verdict: almost
- Key criticisms: (1) MAJOR — the "600-seed paired sweep" numbers were one
  200-seed shard (dict.update last-wins in the summary); (2) MAJOR− —
  "saturates near k≈8–10" unsupported (increments still ~+0.014 at k=9→10);
  (3) eig-only range slip (chr22 is 0.503); (4) "loadings at every scale"
  overreach (no biobank-scale eig-only arm); (5) undisclosed `gauss`
  approximation of the BN frequency draw; (6) unsupported Limitations clause +
  dangling ref (no GPU parser was benchmarked); (7) consolidation defects
  (empty S36 section, wrong S28/S29 labels); (8) Fig 8 caption/panel
  mismatches; (9) s31 3-pop/5-pop not matched; (10) s36 non-member self-excl
  dummy. Verified clean: all 12 ref-basis controls at chance; monotone
  in-sample gap; permutation p<5e-5; density/5-pop ordering; s36 0.486/0.671.

## Round 2 (2026-10-06) — Reviewer backend: pi-review
- Score: 8/10; Verdict: almost
- Both MAJORs genuinely fixed and independently re-verified (independent
  re-pooling reproduced +0.0943 [+0.0925,+0.0962] P=1.0 n=600; eig-only now
  real at biobank scale 0.497–0.511 over 12 cells; gauss disclosure adequate
  and verified immaterial — the 600-seed PN panels use the exact Beta draw).
- Remaining: Fig 8(a) caption false about its own panel (plotted refbasis−AF,
  not a chance-level control); "six shards" miscount; unstated 200-seed cells;
  Limitations "+0.095" vs sec:scale "+0.094"; mixed-column gap endpoints; a
  latent shard-merge inconsistency in the figure path.

## Round 3 (2026-10-06) — Reviewer backend: pi-review — FINAL
- Score: 8.5/10; Verdict: **ready**
- Fixes verified against artifacts + the rendered PNG: Fig 8(a) now plots
  `refbasis−0.5` so the dashed control genuinely sits at ~0 and the caption is
  true; shard/seed counts match the artifact trail (3 shards + one slow run);
  Limitations number reconciled; per-Fst-column endpoints; shard merge made
  first-wins (setdefault) in both summary and figures.
- Residual items are copy-edit only (one endpoint's M-level; "slow cell" vs
  "slow cells"; grid seed count stated indirectly).
- STOP: score 8.5 ≥ 6 AND verdict ∈ {ready, almost} → converged.

### Reviewer Raw Responses
Full verbatim responses for rounds 1–3 are appended in
`review-stage/REVIEWER_MEMORY.md` (reviewer-owned audit trail).

### Reviewer backend
pi-review (model glm-flash, local pi agent), difficulty: medium.
Thread: `a23d84b2238441adb813100537127f22`. run_id: `run_20261006_350bbd97`.
Stop gate: `{"decision":"stop","identity_assurance":"not_required"}`.

### Actions Taken
- `code/s29_negation.py`: per-seed deltas stored; PAIRED_FAST/PAIRED_SLOW split.
- `code/s29_paired_slow.py`, `code/s29_ladder.py`: new (correct pooling pieces).
- `code/gpu_round_summary.py`: `merge_paired` concatenates per-seed deltas;
  eig-only column; S36 schema fixed.
- `code/s27_gpu_scale.py`: eigenvalue-only control added to every cell; shard2
  re-run so all 12 cells carry it.
- `code/s36_phased_ibd.py`: non-member self-exclusion fixed (no `-1` dummy).
- `code/make_figures.py`: fig8(a) plots refbasis−0.5; first-wins shard merge.
- `paper/main.tex`: Section `sec:scale` + Fig 8 caption reworded; Limitations
  gauss/parser disclosures; numbers reconciled.

### Results
- Paired PC-AF gap (pooled): small-$n$ +0.094 [+0.092,+0.096] P=1.0 (600
  seeds); $M{=}8000$ +0.009 P=0.94; $n{=}1600$ −0.014 P=0.21.
- Exact k-ladder monotone, decaying increments; eig-only at chance in all 12
  biobank cells (0.497–0.511), real (0.50–0.53), msprime.
- Paper: `paper/main.pdf`, 14 pp, 0 errors, 0 undefined refs.

### Status
- stopped (converged, positive assessment)
- Difficulty: medium

## Method Description (GPU / scale round, updated)

The leak-measurement harness adds a GPU kernel layer (`code/common_gpu.py`): a
host-resident int8 panel simulator (exact Balding–Nichols `Beta` draw for
small panels; a Gaussian diffusion approximation of the per-population
frequency for biobank-scale cells), streamed AF-LRT scoring, streamed
projection-distance scoring, and a blocked exact-Gram PCA (`torch.linalg.eigh`,
TF32 off). Data flow: host int8 panel -> marker blocks streamed to GPU ->
n×n Gram / attack scores -> AUROC -> `results/raw/*.json` -> tables/figures.
A new experiment family (s27–s36) locates the PC-vs-AF boundary in
(study-size × marker-count × Fst) at biobank dimensions, with a
reference-fitted basis and an eigenvalue-only control in every cell; adds
stratified label-permutation nulls and a per-seed-pooled paired sweep; runs an
exact k=1..10 ladder on sim/real/msprime; sweeps real marker density (chr22/1/2)
across 20× and all five 1000G super-populations; measures the phased
haplotype/IBD channel with planted-IBD and self-exclusion controls; and records
the VCF-parse cost as a CPU-bound negative result. Aggregation is
`code/gpu_round_summary.py` -> `refine-logs/GPU_ROUND_RESULTS.md`; the artifact
map is `refine-logs/PAPER_INDEX.md`. CPU analyses use the project env; the GPU
kernels use the system torch python. All steps are seeded and traceable 1:1.

---

# AUTO_REVIEW — run_20261006_f68b0cd6 (Round 11 — limitation closure)

Backend: pi-review (local glm-flash), difficulty: medium, MAX_ROUNDS=4.
Artifacts under review: `paper/main.tex` / `main.pdf` (15 pp, round-11 revision),
`code/s37_realsynth_scale.py`, `code/s38_privacy_metric.py`, `code/enforce.py`,
`code/s39_enforcement.py`, `code/s40_bundle.py`, `code/s41_phased_ibd_genome.py`,
`code/s42_missingness_realistic.py`, `code/semi_synthetic.py`,
`results/raw/s3[7-9]*.json`, `results/raw/s4[0-2]*.json`,
`refine-logs/EXPERIMENT_RESULTS.md`, `EXPERIMENT_AUDIT.md`.

## Round 1 (2026-10-06) — Reviewer backend: pi-review
- Score: 7.5/10; Verdict: almost. Thread `d3997a6cbc444b8faa219526e48ddd97`.
- Verified all six round-11 artifact groups match the paper to quoted precision;
  stale "remains argued, not built" removed; ~20 round-8/9/10 regression checks
  pass. Weaknesses W1–W8: A*-vs-2AUC−1 direction inverted; eps_lb≈24 clip-driven
  (n_neg=900 not 300); §replication future-work contradiction; Limitations
  cross-ref; saturated base-missingness wording; self-contradictory "full marker
  density"; stale design-doc grid + CV; dead import/docstring.

## Round 2 (2026-10-06) — Reviewer backend: pi-review
- Score: 8.0/10; Verdict: almost.
- W1–W5, W8 fixed; W6 partially; W7 mostly. New residual R1–R5: paper's
  "0.53/4.2 for full-precision AF" was af_2dp (af_full 0.55/4.1); results-doc
  CoV line garbled; tracker S38 direction shorthand stale; s41 still "at full
  marker density"; eps_lb_*_valid=False undisclosed.

## Round 3 (2026-10-06) — Reviewer backend: pi-review — FINAL
- Score: 8.5/10; Verdict: **ready**.
- All W1–W8 and R1–R5 verified fixed against the artifacts; no new defect.
  Regression: paper diff = the reviewed fix blocks only; s37–s42 JSON unaffected;
  15 pp, 0 errors.
- STOP: score 8.5 ≥ 6 AND verdict ∈ {ready, almost} → converged.

### Reviewer Raw Responses
Full verbatim responses for rounds 1–3 are appended in
`review-stage/REVIEWER_MEMORY.md`.

### Reviewer backend
pi-review (model glm-flash, local pi agent), difficulty: medium.
Thread: `d3997a6cbc444b8faa219526e48ddd97`. run_id: `run_20261006_f68b0cd6`.
Stop gate R3: `{"decision":"stop","identity_assurance":"not_required"}`.
Traces: `.aris/traces/auto-review-loop/2026-10-07_run01/`.

### Actions Taken
- `code/s38_privacy_metric.py`: `eps_lower_bound` takes `n_pos, n_neg` and floors
  the false-negative rate at `1/n_pos` (was `1e-12`, clip-driven); dead import
  removed.
- `paper/main.tex`: two-branch eps_lb formula + n_pos=n_neg=900 + ordering-only
  caveat; A*-vs-2AUC−1 direction corrected; §replication IBD sentence redirected
  to §scale; Limitations cohort cross-ref → sec:scale; s41 "uniform genome-wide
  density"; base-missingness wording "stays saturated at ρ=1".
- `code/governance.py`: docstring / TCB note (binding built in `code/enforce.py`).
- `refine-logs/EXPERIMENT_RESULTS.md`, `refine-logs/ROUND11_DESIGN.md`,
  `refine-logs/EXPERIMENT_PLAN.md`: doc corrections (CoV, delivered grid, removed
  MLP/anc-only promises, eps validity note).

### Results
- Round-11 experiments (s37–s42) all DONE; see `refine-logs/EXPERIMENT_RESULTS.md`.
- `paper/main.pdf`, 15 pp, 0 errors, 0 undefined refs.

### Status
- stopped (converged, positive assessment)
- Difficulty: medium

## Method Description (round 11 — limitation closure)

The measurement harness is extended with four governance/measurement additions.
(1) An absolute per-individual privacy metric (`code/s38_privacy_metric.py`)
re-expresses each operating point as the membership advantage A* (=KS), a
per-individual posterior, and an empirical f-DP lower bound with both trade-off
branches and a sample-resolution floor. (2) A built sole-egress binding
(`code/enforce.py`): the study computation runs as a child process whose only
writable channel is one pipe descriptor wired to the parent `EgressMonitor`,
with a PEP-578 audit hook denying write-opens/sockets/subprocess/ctypes; the
bypass is demonstrated against the plaintext prototype then closed
(`code/s39_enforcement.py`). (3) A joint release-bundle attacker
(`code/s40_bundle.py`): out-of-fold logistic over pooled AF, per-site AF, and
the call-rate table, versus the strongest single object, with a structural-zero
noise control and leave-one-out attribution. (4) Two extended second channels:
multi-chromosome physical-length phased IBD (`code/genetic_map.py` +
`code/s41_phased_ibd_genome.py`, cM via the 1000G IMPUTE2 map, planted-IBD
control) and a batch-composition missingness channel on the real 1000G skeleton
(`code/s42_missingness_realistic.py`, monotone in the batch parameter rho with a
rho=0 null). A biobank-scale semi-synthetic cohort with real haplotype structure
is built by Li–Stephens mosaic copying (`code/semi_synthetic.py` +
`code/s37_realsynth_scale.py`) and passed through the s27 attack battery with LD
decay, relatedness, AF-spectrum and eigenvalue-fraction controls. Data flow:
generated or real panels -> attack scores -> AUROC/absolute metrics ->
`results/raw/*.json` -> paper. All steps are seeded and traceable 1:1.

---

# AUTO_REVIEW — run_20261007_3a710cda (Round 12 — residual-limitation closure)

## Round 12 (2026-10-07) — Reviewer backend: pi-review
- Score: 8.5/10; Verdict: almost.
- Verified all four round-12 modules are real: (2) the native-code test genuinely
  raises on `ctypes.dlopen` (`reason='ctypes.dlopen'`, `child_failed=False`, the
  hook is the only `blocked=True` path); (3) s46's uncoarsened comparator uses the
  same full-precision per-site AF as the bundle, and `single_coarse` (0.94707)
  reproduces s40's max_single exactly; (4) s44 real/cohort are real data with
  controls at chance and an honest ordering-only report on real chr22. All
  round-12 paper numbers match JSON keys.
- Weaknesses W12-1..W12-5 (stale "five of six cells" count; results-doc LOO drop
  .076→.233; B4 "S47" tag; missing M=3000 clause; "makes valid" → "admits valid").

## Round 13 (2026-10-07) — Reviewer backend: pi-review — FINAL
- Score: 8.5/10; Verdict: **ready**.
- W12-1..W12-5 all verified landed against the artifacts; no new defect.
  Regression: 16 pp, 0 errors; round-11 and round-12 numbers all preserved.
- STOP: score 8.5 ≥ 6 AND verdict ∈ {ready, almost} → converged.

### Reviewer Raw Responses
Full verbatim responses for rounds 12–13 are appended in
`review-stage/REVIEWER_MEMORY.md`.

### Reviewer backend
pi-review (model glm-flash, local pi agent), difficulty: medium.
Thread: `d3997a6cbc444b8faa219526e48ddd97`. run_id: `run_20261007_3a710cda`.
Stop gate R13: `{"decision":"stop","identity_assurance":"not_required"}`.
Traces: `.aris/traces/auto-review-loop/2026-10-07_run01/004..005`.

### Actions Taken (round 12)
- `code/s38_privacy_metric.py`: `eps_lower_bound` gains `fnr_mode` (resolution /
  rule3); `metrics_from_scores` emits both.
- `code/s44_individual_real.py`: new — s38 estimators on real 1000G chr22 +
  an s37 cohort.
- `code/s46_bundle_manifest2.py`: new — second manifest + uncoarsened single.
- `code/enforce.py` + `code/s39_enforcement.py`: native-code (`ctypes.CDLL`)
  injection; `native_egress_blocked_frac`.
- `code/s37_realsynth_scale.py`: nsp=16000 cell.
- `paper/main.tex`: native-boundary sentence, per-individual real/cohort
  paragraph, second-manifest bundle sentence, nsp=16000 + M=3000, Limitations
  rewritten (16 pp).
- `refine-logs/EXPERIMENT_RESULTS.md`, `EXPERIMENT_TRACKER.md`,
  `refine-logs/EXPERIMENT_PLAN_R12.md` added.

### Results
- s44: real chr22 pc A* .27 < af .35; cohort af .14 > pc .04; controls ≤ .07;
  rule-3 eps valid on cohort, ordering-only on real.
- s46: bundle .973 vs uncoarsened single .964 (+.008, 5/5).
- s39: native-egress blocked 3/3.
- s37: nsp=16000 af .529 vs pc .503 (delta −.026), controls pass.
- `paper/main.pdf`, 16 pp, 0 errors, 0 undefined refs.

### Status
- stopped (converged, positive assessment)
- Difficulty: medium

---

# AUTO_REVIEW — run_20261007_66850dac (Round 13 — real-data structure + size)

## Round 13 (2026-10-07) — Reviewer backend: pi-review
- Score: 8.0/10; Verdict: almost.
- Real-genome ladder (s48) validated: controls real + at chance, s49 held-out-half
  transfer has no leak. But the "measured-structure law" framing over-claimed:
  W13-1 a cited cross-population range hid the chr2 3-pop +0.023 cell; W13-2
  "monotone in structure" contradicted by two rising steps; W13-3 structure
  confounded with study size, uncontrolled; W13-4 the k/n≈1/3 regime undisclosed;
  W13-5 fragmented §boundary sentence, imprecise control ranges, wrong B2 chr.

## Round 15 (2026-10-07) — Reviewer backend: pi-review — FINAL
- Score: 8.5/10; Verdict: **ready**.
- Added the decisive scale-matched rerun (`s48_scalematched.json`, per-group n=20):
  at matched small n the PC channel dominates even across populations (3pop +0.161,
  5pop +0.111), so the AF-first switch needs larger n *and* structure — the
  real-data ordering is governed jointly by study size and structure, the same
  small-(n×M) boundary the simulator locates. Reframed in all five places
  (abstract/Contributions/§scale/Discussion/Limitations); W13-1..W13-5 fixed;
  17 pp, 0 errors.
- STOP: score 8.5 ≥ 6 AND verdict ∈ {ready, almost} → converged.

### Reviewer Raw Responses
Rounds 13/15 verbatim in `review-stage/REVIEWER_MEMORY.md`.

### Reviewer backend
pi-review (glm-flash, local pi agent), difficulty: medium.
Thread: `d3997a6cbc444b8faa219526e48ddd97`. run_id: `run_20261007_66850dac`.
Stop gate R15: `{"decision":"stop","identity_assurance":"not_required"}`.
Traces: `.aris/traces/auto-review-loop/2026-10-07_run01/006..007`.

### Actions Taken (round 13)
- `code/s48_real_structure.py`, `code/s49_real_mechanism.py`,
  `code/s50_real_multichrom.py` (new) + `s48_scalematched.json` (via s48 cohort_cell per=20).
- `paper/main.tex`: abstract, Contributions, §scale ("Where the ordering lives on
  real genomes" + scale-matched paragraph), Discussion, Limitations reframed to
  "study size and structure jointly"; the false "monotone"/"below the real-data
  level" wording removed; k/n disclosed; control ranges corrected.
- `refine-logs/EXPERIMENT_RESULTS.md`, `refine-logs/EXPERIMENT_PLAN_REAL.md`.

### Results
- s48: within-population real gap +0.050…+0.134 (12 cohorts, n_study 28–36);
  2-pop +0.112…+0.172; 3/5-pop −0.005/−0.030 (n 450/715).
- s48 scale-matched (per=20): 3pop +0.161, 5pop +0.111 — PC dominates across
  populations at matched small n.
- s49: top-3 0.88–0.96, full 0.99–1.00, controls 0.55–0.60, oos 0.55–0.58,
  same-half 0.998–1.000.
- s50: chr1/chr2 within +0.051…+0.086; across −0.039…+0.023.
- `paper/main.pdf`, 17 pp, 0 errors.

### Status
- stopped (converged, positive assessment)
- Difficulty: medium

---

# AUTO_REVIEW — run_20261008_0d685cee (Round 17 — limitations re-verify)

## Round 17 (2026-10-08) — Reviewer backend: pi-review — FINAL
- Score: 8.5/10; Verdict: almost.
- Mechanical scan of 122 cited numbers against their JSON keys: 120 pass, 2
  single-character mismatches (W17-1, W17-2). One (W17-2, eigenvalue-only range
  lower bound) was real and fixed (0.55→0.54, s49 min 0.5447); the other (W17-1)
  was a reviewer misread — the paper already reads $0.13$ for the s44 cohort A*
  (JSON 0.1328), so no change.
- All prior-round numbers still verified; paper byte-identical to the round-16
  state apart from this one edit; 17 pp, 0 errors.
- STOP: score 8.5 ≥ 6 AND verdict ∈ {ready, almost} → converged.

### Remaining limitations (all disclosed scope bounds, no blocker)
L1 no real biobank-scale genome cohort (largest external-validity bound; the
real within-population cohorts cap at n≈28–36); L2 no per-split dispersion on the
real round-13 gaps; L3 `_cal` orientation on the real ladder; L4 one M-axis
contrast on real data; L5 resolution-limited f-DP ε_lb; L6 plaintext governance
prototype + declared TCB; L7 constructed missingness channel; L8 3 seeds/cell in
s37/s41/s42.

### Reviewer backend
pi-review (glm-flash), difficulty: medium. Thread `d3997a6cbc444b8faa219526e48ddd97`.
run_id `run_20261008_0d685cee`. Stop gate: `{"decision":"stop"}`.
Traces: `.aris/traces/auto-review-loop/2026-10-08_run01/001`.

### Actions Taken
- `paper/main.tex`: eigenvalue-only control range lower bound 0.55→0.54 (§governance
  per-individual paragraph, matching s49 min 0.5447). No other change.

### Status
- stopped (converged, positive assessment)
- Difficulty: medium

---

## Round 18b (2026-10-08) — cross-family reviewer switch to Codex (gpt-6-astra)

Per user request the loop was re-run with the cross-family reviewer `mcp__codex__codex` (model gpt-6-astra, effort xhigh, read-only), thread `01a11a42-2795-7c53-be84-19b6858e9dcf`. Reviewer backend: codex.

### Score progression (Codex)
- C1: 4/10 — not ready. 7 findings (s40 missingness-object disconnect; s41 estimator bugs; s38 eps_lb not valid; s37 fit cap + AF-only combined; real-data joint-gating overreach + `_cal` orientation; s42 constructed missingness; enforce stderr path) + number mismatches.
- C2: 5/10 — not ready. S41 crashed on the final-marker-mismatch case (out-of-bounds `cm[sidx]`); orientation fix still evaluated the calibration rows; eps/scope language.
- C3: **7/10 — almost**. "No remaining critical implementation defect in the paths checked." Remaining items were a manuscript consistency pass.

### Actions taken (C1 → C3)
- s40/s46: the attacker now consumes the declared, coarsened released call-rate table (`results/raw/s40_bundle.json`, `s46_bundle_manifest2.json`); excess +0.011 / +0.006, 5/5.
- s41: `_cm_run_vs_one` clipped `sidx` to `m-1`; full 3-chromosome re-run (5055 s); JSON regenerated Oct 8 17:37. Real LSR 0.512/0.520/0.471, planted 0.769/0.766/0.749.
- s38: ROC tie-grouped; eps described as an in-sample, resolution-limited summary; binormal-KS reference added then the K-shape inference removed.
- s37: `combined` = joint AF+PC; `n_fit`/`fit_frac` recorded; nsp=12000 re-run to 3 seeds; `code/s37_fitsens.py` (n_fit 5000 vs 20000) shows fit-size-stable gap.
- s48/s49/s50/scalematched: held-out-half orientation, calibration rows excluded, shared cal_mask for paired AF/PC−AF; s49 saves per-field std; all re-run.
- s42 described as a generative sensitivity curve; enforce.py child stderr → /dev/null; s39 re-run all 1.0.
- Paper scope propagated: site-size-dependent crossover; tested-attack scoping for the eigenvalue/IBD negatives; manifest-permitted/coarsened wording; all numbers re-read from the regenerated JSON.

### Status
- Score 7/10, verdict **almost** → satisfies the loop stop condition (score ≥ 6 AND verdict ∈ {ready, almost}).
- Reviewer backend: codex (cross-family, OpenAI). Executor: deepseek-flash (caller-declared). independence: caller_declared / unverified.
- Stopping.

---

## Round 19 (2026-10-08) — Codex (gpt-6-astra), thread 01a11a42-2795-7c53-be84-19b6858e9dcf

### Assessment (Summary)
- Score: 7/10
- Verdict: **ready** (READY for submission: Yes)
- Key criticisms: three nonblocking copyedits only (strict "monotone" residue in abstract/contributions/limitations; one surviving "f-DP bound"; claims-doc "n=450/715" phrasing).

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response</summary>

**7/10. READY for submission: Yes.**

I checked the current source, raw JSON, regenerated Figure 7, and its inclusion in the compiled paper. **No remaining critical scientific or implementation weakness was identified.** No additional experiment is needed for the present, scoped claims.

The corrections are verified:
- **Figure 7:** the old global-crossover annotation is gone. Its three nominal crossing bands match S21 exactly, and the compiled paper contains the updated figure and corrected caption.
- **Numerical summaries:** the eigenvalue-control ranges, S50 population subsets, near-parity/AF-first distinction, S38 approximations, and S37 sensitivity wording now match the underlying results.
- **Statistical interpretation:** S38's docstrings and CLI now describe an in-sample, resolution-limited summary. The historical audits are explicitly superseded, and the claims document correctly identifies the ε counterexample.
- **Earlier substantive conclusions stand:** the corrected real-data evaluation supports within-population PC dominance; the native-route denial is genuine; bundle excess remains **+0.011240/+0.006146**, positive in **5/5 seeds** for both manifests.

Three small wording residues remain, ranked below. These are **nonblocking copyedits**, not critical weaknesses:
1. Strict "monotone" survives in the abstract, contributions, and limitations, although the n=300 curve has two small upward steps. Minimum fix: use "approximately decreasing".
2. "f-DP bound" survives once. Minimum fix: replace it with "trade-off summary".
3. The claims document retains "switch … at n=450/715." Minimum fix: synchronize to "near parity at 450; AF-first at 715."

The score remains 7 because this is a credible, useful measurement contribution with bounded generality: selected attacks, small real-data cohorts, and semi-synthetic scaling. Those limitations are now sufficiently explicit. This final pass verified artifact consistency; it did not repeat the expensive full experiment suite.

</details>

### Actions Taken
- main.tex: "monotone in measured structure" → "approximately decreasing" (abstract, contributions, §limitations).
- main.tex §governance: remaining "the f-DP bound finite" → "the trade-off summary finite".
- CLAIMS_FROM_RESULTS.md: "appears only at n=450/715" → "near parity at n=450 and AF-first at n=715".
- Paper recompiled: 17 pages, 0 errors.

### Results
- Score progression this run (Codex cross-family): C1 4/10 → C2 5/10 → C3 7/10 (almost) → **C4 7/10 (ready)**.

### Status
- **Stopping** — score 7/10 AND verdict "ready" satisfies the stop condition.
- Reviewer backend: codex (cross-family, OpenAI gpt-6-astra). Executor: deepseek-flash (caller-declared). independence: caller_declared / unverified.
- Difficulty: medium.

---

## Final Summary — run_20261008_d9f3a1b2 (Codex cross-family loop)

| Round | Backend | Score | Verdict | Outcome |
|------|---------|------:|---------|---------|
| C1 | codex/gpt-6-astra | 4.0 | not ready | 7 code-level findings (s40 missingness object, s41 estimator, s38 eps_lb, s37 combined/fit cap, real-data gating, s42, enforce stderr) |
| C2 | codex/gpt-6-astra | 5.0 | not ready | s41 crashed on the final-marker-mismatch case; orientation still scored the calibration rows |
| C3 | codex/gpt-6-astra | 7.0 | almost | implementation paths clean; manuscript consistency pass |
| C4 | codex/gpt-6-astra | 7.0 | **ready** | no critical weakness; three copyedits applied |

Final verdict: **ready** (score 7/10 ≥ 6 AND verdict "ready"). Reviewer backend: codex
(cross-family, OpenAI gpt-6-astra). Executor: deepseek-flash (caller-declared).
Independence: caller_declared / unverified (no host-event attestation on this path).

Remaining limitations (disclosed in the paper, not blockers for the scoped claims):
- real-data arm is small per-population cohorts; the PC>AF within-population result is an
  operating-point result, and the size-vs-composition dependence is a qualified
  interpretation (total n grows across the scale-matched ladder);
- the phased haplotype channel is bounded only for the tested longest-shared-run/IBS
  estimators, not in general;
- the eigenvalue-only negative is a random-basis projection ablation;
- the ε trade-off summary is an in-sample, resolution-limited descriptive quantity, not a
  statistically valid f-DP bound.

## Method Description

The paper measures \emph{membership-inference leakage} of the objects a multi-site GWAS
releases, and builds a governance layer over those objects. Given a phased or dosage
genotype panel, a study set is split three ways into disjoint \emph{study members},
attacker \emph{reference}, and independent \emph{non-members}; every attack is a
per-individual score that must rank members above non-members (AUROC). The attack
families are: an allele-frequency likelihood-ratio (pooled, per-site/stratum, at several
release precisions); a projection-distance attack on the study-fitted PCA spectrum
(loadings vs a random-basis / reference-fitted-basis control, and an out-of-sample
transfer check); a haplotype/LD and a phase-IBD longest-shared-run channel measured in
centiMorgans; and a membership-correlated per-SNP call-rate channel. A joint out-of-fold
logistic attacker measures the multi-object \emph{bundle} against its strongest single
object. Results are reported on a parametric Balding--Nichols simulator, an msprime
recombination-realistic demography, real 1000G Phase-3 chromosomes, and a Li--Stephens
mosaic biobank-scale cohort (GPU Gram-PCA).

The governance layer: a declarative \emph{reveal manifest} names each quantity and its
declared precision; an \emph{egress monitor} applies the precision and blocks undeclared
quantities; a hash-chained \emph{ledger} records cross-run refinements. The binding is a
process-isolation runtime (PEP-578 audit hook denying write-opens/sockets/subprocesses/
ctypes, one sanctioned pipe to the monitor, child stdout/stderr discarded). Data flow:
panel $\rightarrow$ three-way split $\rightarrow$ attack scores per release object
$\rightarrow$ AUROC / per-individual $A^\ast,\rho,\varepsilon$ summary $\rightarrow$
manifest-governed release with monitor + ledger.

---

## Round 20 (2026-10-09) — fresh independent Codex thread (gpt-6-astra)

Round backend: codex (cross-family). New thread `01a11b98-9533-72c3-855e-f9f6999ac82e` (independent of the prior thread).

### Assessment (Summary)
- Score: 5/10
- Verdict: not ready
- Key criticisms (all fixed in round C6): (1) novelty vs Zari PSD-2022 PCA MI; (2) missingness result is pre-release info (candidate-only baseline matches; S40 contributor split); (3) s41 releases haplotypes then removes own → direct match AUROC 1.0, relabel; (4) semi-synthetic relatedness control on artificial phase (actual 2.2–2.5× donor); (5) af_by_group dict raises TypeError in the monitor; (6) fig7/s24 λGC-only 0.00 → 0.50; plus: topk keeps largest p-values; ε-ordering statement.

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response</summary>

**Score: 5/10; Verdict: not ready.**

The core finding—membership leakage from study-fitted PC loadings—survives this review. Several supporting claims do not. ...
[See REVIEWER_MEMORY.md round C5 for the full verbatim text.]

</details>

### Actions Taken (round C6)
- Novelty: cited Zari et al. (PSD 2022) and reframed the contribution around the controlled GWAS comparison and the boundary measurements; the release object is named as the loading matrix.
- Missingness: added the candidate-only baseline to s42 (own call count, no released or reference table): it matches the released-table attack curve (0.540/0.664/0.760/0.884/1.000 vs 0.502/0.641/0.762/0.891/1.000), so the paragraph is reframed as batch-associated discrimination, not release exposure. s40 now reports the non-contributor-member M3 AUROC (0.511) alongside the all-member 0.705.
- s41: added `code/s41_directmatch.py` (verbatim self-match AUROC 1.00) and relabeled the leave-self-out LSR arm as a relatedness diagnostic; the raw-haplotype release is stated as trivially leaking.
- Relatedness: `semi_synthetic.generate_cohort(..., return_hap=True)` + `code/s37_phase_control.py`: actual copied-phase relatedness is 2.2–2.5× the donor library; text narrowed to the parsimony-proxy statistic.
- Monitor: `apply_precision`/`_value_hash` handle dict-valued objects; `EgressMonitor.release` now delivers per-group AF at full precision and 2 dp with the ledger chained (s39 re-run: pergroup_full 1.0 / pergroup_2dp 1.0 / chained 1.0).
- Fig7/s24: `s24_crossrun` re-run (λGC-only 0.50 baseline); fig7 regenerated; `topk` gained `order="small"` for the p-value table; the ε-ordering sentence now names the AF full-vs-2dp counterexample.

### Results
- Re-ran s24, s39, s40, s42, s46; new probes s41_directmatch, s37_phase_control. Paper 18 pages, 0 errors, 0 undefined refs.

### Status
- Continuing to the next review round.
- Difficulty: medium.

---

## Round 21 (2026-10-09) — Codex (gpt-6-astra), same fresh thread 01a11b98-9533-72c3-855e-f9f6999ac82e

### Assessment (Summary)
- Score: 7/10
- Verdict: **ready** (no remaining critical correctness defect)
- After C7 fixes, three nonblocking terminology items remained: incomplete loading-matrix terminology; a broader "bound"/"governed reading" wording; a transport sentence.

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response (round 21)</summary>

**Score: 7/10; Verdict: ready.** The substantive C7 fixes hold up; no remaining critical correctness defect. Verified the conditional comparisons (S40 0.951602→0.959446, +0.007844; S46 0.965482→0.971084, +0.005602), actual-phase relatedness disclosure, missingness-as-confounding-control, raw-haplotype vs leave-self-out distinction, per-group AF monitor delivery, Fig7 λGC 0.50 baseline. Remaining nonblocking: (1) loading-matrix terminology incomplete at contributions/table/discussion; (2) haplotype paragraph still says "bound"; "governed reading"; bundle conditional gain overclaims; (3) one transport sentence inaccurate.

</details>

### Actions Taken (round 22)
- Terminology: contributions and discussion now say "study-fitted PC loading matrix"; table row "PC loadings".
- Haplotype paragraph: "a characterization of the tested IBS statistic and the leave-self-out LSR diagnostic, not a bound on the haplotype channel as a whole".
- Removed the "governed reading" phrase; bundle conditional gain scoped to "the AF channels plus the candidate's total call count (a specific comparator, not all candidate-side information)".
- Transport sentence corrected: "does not yet emit per-group AF through the monitored release path (per-group AF is serialized and restored in the parent, but is omitted from the isolated runtime's monitored release loop)".

### Status
- Stopping. Score 7/10 AND verdict "ready" satisfies the stop condition.
- Reviewer backend: codex (gpt-6-astra). Paper: 19 pages, 0 errors, 0 undefined refs.

---

## Round 22 (2026-10-09) — fresh independent Codex thread 01a12007 (gpt-6-astra)

Round backend: codex (cross-family). New thread `01a12007-850b-7ba3-abe6-fb13c35ef44b` (independent of 01a11b98).

### Assessment (Summary)
- Score: 5/10
- Verdict: not ready
- Key criticisms: (1) the central PC-vs-AF comparison supplies calibration side information unequally — a conditioned-AF baseline (balanced-half removal of the known calibration members, self-inclusive residual) beats PC 25/25; (2) the sole-egress binding permitted ordinary Python UDP `sendto`; (3) "largest canonical diagnostic channel" overreaches (M1/M4 are per-release constants, M3 a random vector); (4) the "realistic bundle" mixes contributor populations; (5) s25 QC index bug + a prescribed stopping asymmetry; (6) Table 1 bootstrap label, legacy `max(AUC,1−AUC)` orientation, one-sided-vs-two-sided KS; (7) figure/field misdescriptions (fig5/6 sizes, fig8(a) Fst confound, s27 `combined` AF-only, s41 density wording); (8) two overstrong mechanistic statements (GPU LD construction, eigenvalue-ratio control).

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response</summary>

[See REVIEWER_MEMORY.md round C8 for the full verbatim text.]

</details>

### Actions Taken (round C9)
- Conditioned-AF: reproduced the reviewer's 0.9407 (self-inclusive residual, 25/25) and diagnosed it as the paper's self-reference artifact; added the leave-own-out residual (0.429, below chance, PC wins 25/25) and reported all three residuals in a new paper paragraph. The PC>AF ordering is not an artifact of asymmetric side information.
- UDP: hook denies `socket.__new__`/`sendto`/`sendmsg`/`send`; a UDP probe in `enforce.py`/`s39` reports `udp_egress_blocked_frac`=1.0 (3/3); paper states the closed route.
- Wording: abstract/contributions scoped to the two implemented non-degenerate scores; M1/M4/M3 named as limited controls.
- Bundle: paper calls it a constructed bundle and discloses the contributor sets; the out-of-fold known-membership attacker is stated.
- S25: QC index mapping fixed.
- Metrics: Table 1 caption "mean of the per-seed within-run bootstrap intervals"; s21 orientation named direction-agnostic with the null bias disclosed; s38 `membership_advantage` is the two-sided KS.
- Figures: fig5/fig6 axis labels; fig8(a) one curve per (size, Fst); s27 `combined` joint AF+PC; s41 "uniformly strided".
- Mechanism: GPU LD-construction difference disclosed; eigenvalue-ratio caveat added.

### Results
- Re-ran s10 (conditioned/self-ref/self-ref-loo, 25 seeds), s38, s39, s24; regenerated figures; paper 20 pages, 0 errors, 0 undefined refs.

### Status
- Continuing to round C10 on the same thread.
- Difficulty: medium.

---

## Round 23 (2026-10-09) — Codex (gpt-6-astra), same thread 01a12007

### Assessment (Summary)
- Score: 5/10
- Verdict: not ready
- The leave-own-out result was rejected as label-dependent; the reviewer reproduced a class-symmetric candidate subtraction (0.940966) that also beats PC. Other C8 items closed except a factual error about which individuals remain in the residual, S27 artifact refresh, and two label/size issues.

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response</summary>

[See REVIEWER_MEMORY.md round C9 for the full verbatim text.]

</details>

### Actions Taken (round C10)
- Accepted the reviewer's central point: the balanced known-member residual (`af_conditioned` = 0.941) is a valid conditioned AF attack, not a self-reference error. The factual error (the residual removes the *train* members, keeping every evaluated member) is corrected.
- Added the symmetric PC arm `pc_conditioned`: PCA refit on the same remaining panel, one predefined orientation. Mean 0.9636 ± 0.006, above the conditioned AF on 25/25 seeds (paired +0.023, 95% CI [+0.019,+0.028]). This becomes the ordering claim.
- Renamed the member-only leave-own-out channel `af_loo_oracle`, labeled an oracle ablation.
- S27 `combined` re-run underway (joint AF+PC); the previous AF-only values were the artifact.
- fig6(d) axis "study individuals"; fig5(c) ticks show study individuals; Table 1 caption "TPR@1% bounds".
- Abstract eigenvalue sentence, discussion "constructed bundle", bundle attacker split, UDP probe description propagated.

### Status
- Continuing to round C11 on the same thread.
- Difficulty: medium.

---

## Round 24 (2026-10-09) — Codex (gpt-6-astra), same thread 01a12007

### Assessment (Summary)
- Score: 5/10
- Verdict: not ready
- The reviewer accepted the balanced conditioned-AF interpretation but rejected the framing that `pc_conditioned` (a PCA refit on the remaining members' genotypes) restores PC superiority under a matched attacker-information model. Two minor artifact/label items also remained.

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response</summary>

[See REVIEWER_MEMORY.md round C10 for the full verbatim text.]

</details>

### Actions Taken (round C11)
- Paper now states the as-released ordering honestly: conditioned AF (0.941) is the strongest implemented attack against the objects as released, above the full-study PC (0.879) and the AF+PC combination (0.935); the headline PC>AF increment is against the unconditioned AF baseline. `pc_conditioned` (0.964) is explicitly framed as a smaller-cohort curator-refitted release comparison, not an attack on the original release. Abstract, contributions, CLAIMS updated.
- Removed the misleading auxiliary S27 files; regenerated `gpu_round_summary.json` from the corrected shards.
- Fig5 counts use `3*int(0.3*n)`; the `pc_conditioned` orientation is described as reused from training calibration.

### Status
- Continuing to round C12 on the same thread.
- Difficulty: medium.

---

## Round 25 (2026-10-09) — Codex (gpt-6-astra), same thread 01a12007

### Assessment (Summary)
- Score: 7/10
- Verdict: **ready** (no remaining submission-blocking defect; the comparison scope is accurate)
- The reviewer recalculated the S10 summaries and confirmed the as-released ordering, the curator-refit separation, and the S27/Fig5 fixes. Only two stale inline code comments remained.

### Reviewer Raw Response

<details>
<summary>Click to expand full reviewer response (round 25)</summary>

**Score: 7/10; Verdict: ready.** The manuscript's comparison scope is now accurate. The previous substantive findings are closed. Conditioned AF 0.940730 (original PC 0.878575, 25/25); AF+PC combination 0.935391 (conditioned AF beats it 19/25, a higher mean); curator-refitted PCA 0.963559, correctly separated. Remaining: two stale inline comments in s10_paired.py. "I found no remaining submission-blocking scientific or artifact defect … I would submit it."

</details>

### Actions Taken (final)
- s10_paired.py: the two stale comments replaced ("curator-refitted smaller-cohort release", "training-calibrated sign").

### Status
- Stopping. Score 7/10 AND verdict "ready" satisfies the stop condition.
- Reviewer backend: codex (cross-family, OpenAI gpt-6-astra). Executor: deepseek-flash (caller-declared). Independence: caller_declared / unverified.
- Difficulty: medium.

---

## Final Summary — run_20261009_27c9f1a4 (Codex cross-family loop, thread 01a12007)

| Round | Backend | Score | Verdict | Outcome |
|------|---------|------:|---------|---------|
| C8 (22) | codex/gpt-6-astra | 5.0 | not ready | unequal calibration info in the PC-vs-AF comparison; UDP egress; largest-channel overreach; bundle contributors; s25 index; labels; figures; mechanism |
| C9 (23) | codex/gpt-6-astra | 5.0 | not ready | leave-own-out rejected as label-dependent; balanced conditioned AF reproduced 25/25; S27 artifacts; fig/table labels |
| C10 (24) | codex/gpt-6-astra | 5.0 | not ready | the PCA refit is a curator action, not an attack on the original release; S27 auxiliary files; fig5 counts |
| C11 (25) | codex/gpt-6-astra | 7.0 | almost | scope propagation into summary claims; s10 comments |
| C12 (25) | codex/gpt-6-astra | 7.0 | **ready** | comparison scope accurate; no blocking defect; routine comment cleanup |

Final verdict: **ready** (score 7/10 >= 6 AND verdict "ready"). Reviewer backend: codex
(cross-family, OpenAI gpt-6-astra). Executor: deepseek-flash (caller-declared).
Independence: caller_declared / unverified.

Headline resolution of this loop: the central PC-vs-AF comparison now separates attacker
post-processing from curator refitting. The conditioned AF residual (subtract known
calibration members from the released summary) reaches 0.941, above the full-study PC
(0.879) on 25/25 seeds and above the AF+PC combination (0.935) on 19/25; the headline
PC>AF increment (+0.046; combined +0.103) is against the unconditioned AF baseline. The
PCA refit on the same remaining panel (0.964) is a smaller-cohort curator-refitted release
comparison, not an attack on the original release. The leave-own-out residual (0.429) is a
label-dependent oracle ablation.

Remaining limitations (disclosed in the paper, not blockers for the scoped claims):
- real-data arm is small per-population cohorts; the PC>AF within-population result is an
  operating-point result;
- the phased haplotype channel is bounded only for the tested estimators;
- the eigenvalue-only negative is a random-basis projection ablation;
- the epsilon trade-off summary is descriptive, not a valid f-DP bound.

---

## Round R15-0..R15-3 (2026-10-09/10) — round-15 limitation-closure loop

Backends: Codex/gpt-6-astra (R15-0 thread 01a1217c, score 6/10) and pi-review/glm-flash (R15-1/R15-2, thread fa6ff704). Independent reviewer for the loop = pi-review/glm-flash (cross-family).

| Round | Reviewer | Score | Verdict | Outcome |
|------|---------|------:|---------|---------|
| R15-0 | codex/gpt-6-astra | 6.0 | almost | L2 validation circular; L1 "pure total-size effect" overstatement; pooled "zero structure" manufactured by labeling; S44 stale after the helper change; minor number fixes |
| R15-1 | pi-review/glm-flash | 7.5 | almost | all four R15-0 findings fixed; W1-W6 consistency items (s53 not in the body; stale docs; claims lag; wording; comment; overfull) |
| R15-2 | pi-review/glm-flash | 8.0 | **ready** | W1-W5 fixed; N1/N2 copyedit nits; no blocking defect |

Round-15 experiments (all real, number-traceable): s51 (real low-structure / fixed-superpopulation-proportion ladder), s52 (two further leave-self-out haplotype statistics + fixed planted control), s53 (fixed + eigenvalue-weighted bases at chance), s38 (floor-convention consistency, inversion example), s37_fitsens (fit-fraction sweep).

### Reviewer Raw Response (R15-2, ready)

<details>
<summary>Click to expand</summary>

**Score: 8/10; Verdict: READY.** Every position checked now states what its artifacts show, including the body positions that argue the punchline. The residual defects are copyedit-class; none touches a claim, number, or experiment.

</details>

### Final state
- Verdict **ready** (8/10 >= 6 AND verdict "ready"). Reviewer backend: pi-review/glm-flash (cross-family, different family from the deepseek-flash executor).
- Paper: 21 pages, 0 errors, 0 undefined refs.
- EXPERIMENT_AUDIT (round 15): PASS.
- Residual non-blocking: two overfull hboxes; the disclosed scope bounds (zero-structure real large-n unobtainable inside 1000G; semi-synthetic scale arm; tested-attack scoping).

---

## Extension loop (2026-10-10) — chromsome extension of the real-data arms

Backend: pi-review/glm-flash (cross-family). Extended s51 (real size ladder) to chr22+chr1+chr2, s52 (second haplotype statistics) to chr22+chr1+chr2, s53 (fixed eigenvalue-informed bases) to simulated + chr22/chr1/chr2 × EUR/EAS/AFR.

| Round | Score | Verdict | Outcome |
|------|------:|---------|---------|
| EXT-1 | 8.5 | ready | chr2 wording/bad-control checks pass; W1-W4 (tracker wording; seeds paragraph; rung noise; s53 coverage) |
| EXT-2 | 9.0 | **ready** | W1-W4 fixed; s53 extension clean; N1 (range attribution) and N2 (docstring) — copyedit |

### Final state
- Verdict **ready** (9/10 >= 6 AND verdict "ready"). Reviewer: pi-review/glm-flash (cross-family).
- Paper: 21 pages, 0 errors, 0 undefined refs.
- Extension results: s51 pooled chr2 PC>AF +0.111..+0.181, ladder near-parity at the top rung; s52 planted control .772-.798 on chr2; s53 fixed bases 0.464-0.530 across 9 real cells vs full loadings 0.79-0.90.

---

## Round 18-1 (2026-10-14) — run_20261014_a3e66379

Backend: pi-review / glm-flash (cross-family; thread 2da7cc2b7a5e4840b085d168354c5077).

### Assessment (Summary)
- Score: 7/10
- Verdict: almost (ready after text fixes; no new experiments)
- Key criticisms: W1 s42 candidate-only wording (reversed at rho=0.75; violated at rho=0), W2 s53 range min, W3 s60 tautological clause, W4 Fig 8(a)/(b) attribution, W5 s55 PC saturation undisclosed, W6 stale tracker rows, W7 minor ε range + overfull.

### Actions Taken
- W1–W7 all fixed in paper/main.tex + EXPERIMENT_TRACKER.md (see REVIEWER_MEMORY.md).
- Recompiled: 22 pages, 0 errors.

### Status
- Continuing to round 2.
- Difficulty: medium.

## Round 18-2 (2026-10-14)

### Assessment
- Score: 7.5/10
- Verdict: almost
- W1–W7 verified fixed against source JSONs; one new defect N1: abstract claimed the s56 fixed-reference non-degenerate gap was "+0.06 to +0.10" but the AFR full-reference cells span +0.095–0.188.

### Actions Taken
- N1 fixed: abstract now "+0.10 to +0.19"; recompiled 22 pp, 0 errors.

## Round 18-3 (2026-10-14) — FINAL

### Assessment
- Score: 8/10
- Verdict: **ready**
- Reviewer verified N1 fixed against s56 AFR span (+0.0951–0.1876 → "+0.10 to +0.19"); no remaining critical weaknesses.
- STOP CONDITION: score 8 >= 6 AND verdict "ready".

### Final state
- Verdict **ready** (8/10). Reviewer: pi-review/glm-flash (cross-family, thread 2da7cc2b7a5e4840b085d168354c5077). run_id run_20261014_a3e66379.
- Paper: 22 pages, 0 errors, 0 undefined refs.
- EXPERIMENT_AUDIT (round 17): WARN (integrity pass; all WARN items fixed).
- Residual non-blocking: s42 "within 0.021" vs 0.0212 (cosmetic), 5 overfull hboxes, long abstract.

## Method Description
The paper measures the membership leakage of a multi-site GWAS run's released diagnostics and governs them. The measurement uses a controlled MI protocol: a 3-way disjoint split of a real/simulated panel into study members, an attacker reference, and independent non-members; an AF-LRT baseline and a study-fitted PC-projection-distance channel; held-out-half orientation (sign on a calibration half, AUROC on the disjoint other half); and two structure controls (reference-fitted basis, eigenvalue-only random basis). The headline finding is that the released study-fitted PC loading matrix is the stronger intermediate leak in the small-(n×M) regime (0.88 vs 0.84 AUROC), that the leak is the loadings' in-sample-ness (eigenvalue-only at chance, refit-on-fresh-panel at chance, held-out transfer at chance), and that per-site AF is the largest single marginal-frequency channel. The ordering is a boundary phenomenon: PC dominates within a population on real 1000G and across populations at small n, switching to AF-first once study size and composition grow; a GPU-ported harness locates the crossover at biobank scale. The governance layer is a declarative reveal manifest, a process-isolated egress monitor (PEP-578 audit hook, single sanctioned pipe), and a hash-chained cross-run ledger, stress-tested against careless-author bypasses and a native-code boundary. Data flow: cached 1000G panels / simulated panels / GPU-generated cohorts → attack scores → per-channel AUROC/gap tables (results/raw/*.json) → figures and paper claims, with every quoted number pinned to a JSON key.
