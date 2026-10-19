# Experiment Plan (Round 15 — limitation closure for the current 20-page draft)

**Problem**: The round-14 Codex loop converged at 7/10 *ready*, but five
limitations remain explicit in the paper and the claims file. Each is a
measurement-bound gap, not a mechanism gap:
(a) the real-data arm never reaches a **large-n, low-structure** operating point,
so the AF-first switch is only shown on the simulator and on a *semi-synthetic*
cohort;
(b) the **eps metric** is not a valid f-DP bound (`valid=False` on the top
channels) and is not order-preserving with AUROC;
(c) the **phased-haplotype channel** is bounded only for the tested IBS/longest-
shared-run estimators, not in general;
(d) the **eigenvalue-only negative** is a random-basis ablation, not an argument
against any eigenvalue-dependent attack;
(e) the **large-cohort cells** fit the PCA on a subsampled study (<=20000 rows).

**Method Thesis**: Each limitation is closed by re-running an existing estimator
at a new operating point or against a second construction, without new
machinery.
**Date**: 2026-10-09

## Claim Map
| Claim | Why It Matters | Minimum Convincing Evidence | Linked Blocks |
|-------|-----------------|-----------------------------|---------------|
| L1 The AF-first switch is reachable on **real** 1000G at large n in a low-structure cohort | The current large-n cells are simulated or semi-synthetic; a privacy reviewer reads "the real switch is extrapolated" | A real single-super-population cohort at the largest n the 1000G panel supports, gap <= 0, with both basis controls at chance, and a fixed-composition n-ladder holding the group count fixed | B1 |
| L2 The eps summary is not floor-driven and its non-order-preservation is bounded | `valid=False` on the high-AUROC channels leaves the absolute number unquotable | A rule-of-three confidence FNR floor so `valid=True` at some delta on the top channel, plus a stated AUROC-vs-eps inversion example | B2 |
| L3 The haplotype-channel negative is not one estimator | "bounded only for the tested estimators" is an open claim | A second, independent haplotype/LD attack (per-window LD-correlation or identical-by-state share) also fails to separate unrelated members | B3 |
| L4 The eigenvalue-only result is scoped correctly, or a stronger eigenvalue-dependent attack is probed | The random-basis ablation is weak evidence | A **fixed** (non-random) eigenvalue-informed basis is at chance, or an explicit test showing eigenvalues alone do not carry the channel | B4 |
| L5 The subsampled-fit cells match the full-fit result | "subsampled PCA" is a disclosed but unverified approximation | `s37_fitsens` extended to the largest delivered cell, or a fit-fraction sweep at two large cells showing the gap is stable | B5 |

## Paper Storyline
- **Main paper must prove**: L1 (the real large-n low-structure cell). This is
  the strongest remaining upgrade because it turns the cohort arm from
  "semi-synthetic" into a real-genome result on both sides of the switch.
- **Appendix can support**: L2 (eps correction), L3 (second haplotype estimator),
  L4 (eigenvalue scoping), L5 (fit-fraction stability).
- **Experiments intentionally cut**: a UKB-scale real cohort (not obtainable);
  rare-variant/CNV channels (out of scope); a full 22-autosome sweep (chr1/2/22
  suffice).

## Experiment Blocks

### Block 1 (B1): Real large-n, low-structure ladder
- Claim tested: L1.
- Why this block exists: the AF-first switch is currently simulated-only at large
  n. Two real probes (chr22, cached, `scratch/probe_lowncell.py` +
  `scratch/probe_ladder.py`, 2026-10-09) show what is reachable:
  - single super-population (EUR/EAS, measured structure 0.0) at per-group
    n=60/100/150 stays **PC>AF** (gap +0.154/+0.100/+0.086), so 1000G's
    per-population size (~500 samples, ~150-250 per split) is **below** the real
    switch;
  - a **fixed-composition** 3-super-population cohort (measured structure 0.088,
    every rung) crosses at fixed structure: per-group n=40/80/120/160 gives
    gap **+0.120 / +0.042 / +0.007 / −0.026** (total n 120→480), monotone
    decreasing and crossing zero at ntot≈400.
  So the real switch is a total-size effect at fixed composition, and inside
  1000G the low-structure arm cannot reach it.
- Dataset / split / task: real 1000G phase-3, existing `s48.load_for` /
  `cohort_cell`. Arm A: single super-population (EUR, EAS, AFR) at the largest
  per-group n the split allows (report the cap as an honest bound). Arm B: the
  fixed 3-super-population composition ladder per-group n in {40, 80, 120, 160},
  chr22 then chr1/2. 5 split seeds each.
- Compared systems: AF-LRT vs study-fitted PC projection distance, plus the
  reference-fitted-basis and eigenvalue-only controls.
- Metrics: `(af_auc, pc_auc_calibrated, gap, fst_like, refbasis, eigonly,
  n_study, M)`, averaged over split seeds with std.
- Setup: CPU, project env; chr22 cached (2504 x 3391); chr1/chr2 from cache.
- Success criterion: the fixed-composition ladder shows a monotone gap decrease
  that crosses zero as total n grows (already observed on chr22), with both
  controls at chance in every rung; the single-super-population arm reports the
  panel-size bound.
- Failure interpretation: if the fixed-composition ladder does not cross on
  chr1/chr2, the chr22 crossing is chromosome-specific and C6 stays as is.
- Table / figure target: fills the "no real low-structure large-n cell" row with
  the honest bound plus the fixed-composition crossing; a point added to the
  gap-vs-n real panel.
- Priority: MUST-RUN.

### Block 2 (B2): eps metric — confidence floor + non-order example
- Claim tested: L2.
- Why this block exists: `valid=False` on the four high-AUROC channels.
- Dataset / split / task: reuse the s38 pipeline (`code/s38_privacy_metric.py`);
  add a rule-of-three FNR floor `max(1/n_pos, 3/n_pos)` and recompute `valid`.
- Compared systems: current `1/n_pos` floor vs rule-of-three floor, on the
  simulated panel, real chr22, and the cohort.
- Metrics: `eps_lb`, `valid`, `delta`, and the AUROC-vs-eps inversion pair.
- Setup: CPU; a few seconds per cell.
- Success criterion: `valid=True` at some delta on at least the top channel, and
  the inversion example recorded with its two AUROCs and eps values.
- Failure interpretation: if even the rule-of-three floor leaves `valid=False`,
  report the metric as purely descriptive (already the current wording) and drop
  the bound claim entirely.
- Table / figure target: methodology note; no figure.
- Priority: NICE-TO-HAVE (cheap).

### Block 3 (B3): second haplotype/LD estimator
- Claim tested: L3.
- Why this block exists: the haplotype negative currently rests on the
  longest-shared-run / IBS estimators.
- Dataset / split / task: real 1000G phased chr22/chr1; a per-window LD
  correlation (r^2 between the target and the study's tag sets) or an
  identical-by-state allele-share statistic, oriented the same way.
- Compared systems: the existing LSR/IBS statistic vs the new estimator, both on
  the same members/non-members.
- Metrics: per-estimator AUROC with disjoint negatives and the leave-self-out
  diagnostic.
- Setup: CPU; phased caches exist (`real_phased_chr*`).
- Success criterion: the new estimator also fails to separate unrelated members
  (AUROC near chance), with the planted-positive control still recovering.
- Failure interpretation: if the new estimator separates members, the paper's
  "tested estimators only" scoping must be tightened further to name that
  estimator as a positive counterexample.
- Table / figure target: appendix row.
- Priority: NICE-TO-HAVE.

### Block 4 (B4): eigenvalue-dependent attack probe
- Claim tested: L4.
- Why this block exists: "random-basis ablation" is weak evidence against
  eigenvalue-dependent attacks.
- Dataset / split / task: simulated + real within-population; construct a *fixed*
  eigenvalue-informed basis (e.g. the top eigenvectors of the study's covariance
  reconstructed from the released eigenvalues and an assumed isotropic
  off-diagonal), and score with it.
- Compared systems: full loadings vs random-basis vs fixed eigenvalue-informed
  basis.
- Metrics: AUROC.
- Setup: CPU.
- Success criterion: the fixed eigenvalue-informed basis is at chance, matching
  the random-basis ablation, so the negative is not an artifact of randomness.
- Failure interpretation: if the fixed basis separates members, the paper's
  "eigenvalues do not carry the channel" must be strengthened to acknowledge a
  reconstructible eigenvalue-dependent attack.
- Table / figure target: the eigenvalue paragraph.
- Priority: NICE-TO-HAVE.

### Block 5 (B5): fit-fraction stability at the largest cells
- Claim tested: L5.
- Why this block exists: the large-cohort cells subsample the PCA fit.
- Dataset / split / task: reuse `code/s37_fitsens.py` at the largest delivered
  cell(s); sweep `n_fit` in {5000, 10000, 20000, full}.
- Compared systems: the AF-first gap at each fit fraction.
- Metrics: `gap`, `pc_auc`, `af_auc`.
- Setup: GPU (`/opt/miniforge/bin/python3`), one free GPU; the Gram is
  `n_fit x n_fit`.
- Success criterion: the gap changes by less than the per-cell std across fit
  fractions.
- Failure interpretation: if the gap is fit-fraction sensitive, report the
  direction and magnitude and scope the large-cell claim to `n_fit`.
- Table / figure target: the subsampling sentence.
- Priority: NICE-TO-HAVE.

## Run Order and Milestones
| Milestone | Goal | Runs | Decision Gate | Cost | Risk |
|-----------|------|------|---------------|------|------|
| M0 | sanity | re-run B1 chr22 probe (already run) | the real single-pop cells return a value | minutes | chr1/2 parse if not cached |
| M1 | main upgrade | B1 full (EUR/EAS/AFR + fixed-composition ladder, 5 seeds, chr22 then chr1/2) | gap <= 0 at large n, or monotone decrease | CPU-hours | the panel caps single-pop n below the switch |
| M2 | decision | B2 + B5 (cheap) | `valid=True` at some delta; gap stable in n_fit | minutes (B2) + one GPU job (B5) | GPU contention |
| M3 | polish | B3 + B4 | both negatives replicate | CPU-hours | phased-parse time |

## Compute and Data Budget
- Total estimated GPU-hours: ~0.5 (B5 only).
- CPU: B1 is minutes on chr22 (cached), up to ~10 min/chromosome for chr1/2 from
  cache; B3 needs the phased caches (present).
- Data preparation needs: none (all cached).
- Human evaluation needs: none.
- Biggest bottleneck: the 1000G single-super-population size cap (EUR 503 / EAS
  504 / AFR 661), which limits how large the real low-structure cell can be.

## Risks and Mitigations
- The real single-pop large-n cells may stay PC>AF (the pool caps at ~500, and a
  two-way split gives per-group ~250). Mitigation: report the honest bound (the
  switch is above 1000G's per-population size) and lean on the fixed-composition
  ladder to isolate size from composition.
- GPU contention on the shared 8x4090 box. Mitigation: pick the least-loaded GPU,
  cap BLAS threads, run B5 alone.
- B3 phased parse time. Mitigation: use the existing `real_phased_chr*` caches.

## Final Checklist
- [x] Main paper tables covered (B1 adds the missing real low-structure large-n row)
- [x] Novelty isolated (no new mechanism; operating-point closure)
- [x] Simplicity defended (no new machinery)
- [x] Frontier contribution justified or explicitly not claimed (no frontier primitive)
- [x] Nice-to-have runs separated from must-run runs (B1 must-run; B2-B5 nice-to-have)
