# Experiment Plan (Round 15b — remaining disclosed-limitation closure)

**Problem**: Round 15b closed the real-gap uncertainty (s54), the real density
axis (s55) and raised s37/s42 to 5 seeds. The paper (21 pp, 0 errors) still
discloses a small set of scope bounds. Each is a measurement reach gap, not a
mechanism gap. The open items, taken verbatim from `paper/main.tex` Limitations
and `CLAIMS_FROM_RESULTS.md` "What the results do not support":

(a) **No real low-structure large-n cell** — the AF-first switch at $F_{ST}\approx0$
is undetermined on real data; the pooled structure-zero arm still has PC>AF at
$n_{\text{study}}{=}167$–$220$ (1000G pool cap).
(b) **eps_lb is not a statistically valid f-DP bound** and is not order-preserving
with AUROC (AF-full AUROC 0.845 > AF@2dp 0.840 but $\varepsilon$ 4.120 < 4.190).
(c) **The phased-haplotype negative** holds only for the tested leave-self-out / IBS
estimators, not in general. It is now a raw-self-match upper bound plus two
further statistics.
(d) **The eigenvalue-only negative** rests on fixed/random data-independent bases,
not a proof that no eigenvalue-dependent attack exists.
(e) **S37 large cells use a subsampled PCA fit** (≤20000 rows; fit-fraction sweep
now 5-seed-stable).
(f) **Governance is a plaintext, release-object-level prototype** under a stated
threat model (no MPC/SPU; local user with process control outside the model).
(g) **The missingness channel is constructed-by-design** (ρ=0 null at 0.490; the
candidate-only baseline matches).

**Method Thesis**: Every remaining item is closed by (i) re-running an existing
estimator at a new operating point, (ii) adding a second branch of an existing
claim, or (iii) writing a stated scope bound where the data genuinely does not
exist — no new machinery, no new privacy primitive.
**Date**: 2026-10-13

## Claim Map
| Claim | Why It Matters | Minimum Convincing Evidence | Linked Blocks |
|-------|-----------------|-----------------------------|---------------|
| (a) The real AF-first switch is characterised, not merely extrapolated | "the real switch is only semi-synthetic" is the single largest external-validity read | Either a real within-super-population cohort near the switch (best-effort bound) or an honest panel-cap proof that 1000G cannot supply it, with the semi-synthetic arm carrying the large-$n$ side | B1 |
| (b) The ε summary is not floor-driven and its non-order-preservation is bounded | `valid=False` on the top channels leaves the absolute number unquotable | A stated floor convention applied consistently, plus one documented AUROC-vs-ε inversion | B2 (done in r15; only the wording generalisation remains) |
| (c) The haplotype negative is not one estimator | "bounded only for the tested estimators" is open | A second independent statistic (per-window $r^2$ / IBS share) also fails to separate unrelated members, with the raw self-match kept as the total-leak upper bound | B3 |
| (d) The eigenvalue-only result is scoped correctly | The random-basis ablation is weak evidence | A fixed eigenvalue-informed basis at chance, plus an explicit negative statement for the un-probed family | B4 |
| (e) The subsampled-fit cells match the full-fit result | "subsampled PCA" is disclosed but was unverified | A fit-fraction sweep at two large cells showing the gap is stable | B5 |
| (f) A multi-site transport prototype, not an MPC implementation | Governance claim scope | A stated TCB + a process-isolation prototype with a commented sample-protocol boundary | B6 |

## Paper Storyline
- **Main paper must prove**: (a) — the real switch. It is the last row where a
  reviewer can say "the real-data arm never reaches the advertised operating
  point". Everything else is disclosed scope.
- **Appendix can support**: (c) second haplotype estimator; (d) eigenvalue
  scoping; (e) fit-fraction; (b) the ε convention note; (f) the TCB note.
- **Experiments intentionally cut**: a UKB-scale real cohort (not obtainable);
  rare-variant/CNV channels (out of scope); a full 22-autosome sweep (chr1/2/22
  suffice); any f-DP reformulation (that is a different paper).

## Experiment Blocks

### Block 1 (B1): Real switch — pooled structure-zero arm at the true pool maximum + panel-cap bound
- Claim tested: (a).
- Why this block exists: the pooled arm in s51 currently caps at
  $n_{\text{study}}{=}167$–$220$ because it splits the pool three ways. That is
  not the largest real low-structure cohort 1000G supports. A study/reference
  split that is STUDY-HEAVY (e.g. study 60\% / reference 25\% / non-member 15\% of
  the pool) raises the study $n$ to ~300–400 while keeping measured structure
  exactly zero (one group). This is the largest real, structure-zero cohort the
  1000G panel can supply; if it is still PC>AF, the honest statement is "the real
  structure-zero switch is above the 1000G single-population pool".
- Dataset / split / task: real 1000G phase-3 chr22+chr1+chr2, pooled
  single-super-population (EUR/EAS/AFR), study-heavy 3-way split. 5 split seeds.
- Compared systems: AF-LRT vs study-fitted PC projection distance, plus the
  reference-fitted-basis and eigenvalue-only controls.
- Metrics: `af_auc, pc_auc, gap, gap_std, fst_like, refbasis, eigonly, n_study, M`.
- Setup: CPU, project env; cached panels.
- Success criterion: report the largest real structure-zero cohort and its gap;
  either PC>AF (gap published as the panel-cap bound) or a crossover (the real
  switch is exhibited).
- Failure interpretation: if the pooled arm stays PC>AF at every split fraction,
  the real structure-zero switch is unreachable inside 1000G — publish that as a
  measured bound.
- Table / figure target: the "no real low-structure large-n cell" row.
- Priority: MUST-RUN.

### Block 2 (B2): ε convention — global floor statement (analysis only)
- Claim tested: (b).
- Why this block exists: the current wording is per-channel; a reviewer will ask
  for one convention statement. No new experiment — only one sentence in the
  method note that the floor applies on every channel and is not a guarantee.
- Success criterion: the paper has a single ε convention paragraph; no `valid`
  flag is used to separate channels.
- Priority: NICE-TO-HAVE (text only).

### Block 3 (B3): Second haplotype/LD estimator (per-window LD correlation)
- Claim tested: (c).
- Why this block exists: the haplotype negative should not depend on the IBS
  family. Add a per-window LD-correlation statistic (correlation of the target's
  window-wise $r^2$ profile against the release's) to complement the s52
  IBS-share statistics.
- Dataset / split / task: real 1000G phased chr22/chr1/chr2; a planted-IBD
  positive control.
- Success criterion: the new estimator is at chance on the real unrelated arm and
  detects the planted control.
- Failure interpretation: if it separates members, name it as a positive
  counterexample and tighten the scoping.
- Priority: NICE-TO-HAVE.

### Block 4 (B4): Eigenvalue-dependent attack probe (in-simulation only)
- Claim tested: (d).
- Why this block exists: the negative currently rests on data-independent bases.
  Probe the strongest plausible eigenvalue-only attack: reconstruct the study
  covariance from the released eigenvalues under an assumed eigenvector prior
  (e.g. a Toeplitz/LD prior) and project with the reconstructed basis. Fully
  deterministic, in-simulation.
- Success criterion: the reconstructed basis is at chance, or a counterexample is
  found and reported.
- Priority: NICE-TO-HAVE.

### Block 5 (B5): Fit-fraction stability at the largest cells
- Claim tested: (e).
- Status: closed in round 15b (s37_fitsens 5-seed-stable, −0.034 at every
  $n_{\text{fit}}$, full fit included).
- Priority: DONE.

### Block 6 (B6): Governance TCB + transport note
- Claim tested: (f).
- Why this block exists: the governance layer is a plaintext prototype. State
  the TCB and the sample-protocol boundary explicitly; the sole-egress binding is
  realized as process isolation (PEP-578 hook), not MPC.
- Success criterion: a stated TCB paragraph; the native-code boundary tested
  (already 3/3 ctypes denial).
- Priority: NICE-TO-HAVE (text only).

## Run Order and Milestones
| Milestone | Goal | Runs | Decision Gate | Cost | Risk |
|-----------|------|------|---------------|------|------|
| M0 | sanity | B1 pool-fraction probe on chr22 (already cached) | returns a gap | minutes | none |
| M1 | main upgrade | B1 full (EUR/EAS/AFR, study-heavy split, chr22 then chr1/2, 5 seeds) | a crossover, or a stated panel-cap bound | CPU-minutes/cell | the pooled arm may stay PC>AF |
| M2 | decision | B2 + B6 (text only) | convention stated | minutes | none |
| M3 | polish | B3 + B4 | both negatives replicate | CPU-hours (phased parse) | phased-parse time |

## Compute and Data Budget
- Total estimated GPU-hours: 0 (B1/B3/B4 are CPU on cached panels).
- CPU: B1 minutes on cached chr22, up to ~10 min/chromosome for chr1/2 from cache.
- Data preparation needs: none (all cached; `real_phased_*` for B3).
- Human evaluation needs: none.
- Biggest bottleneck: the 1000G single-super-population pool cap (EUR 503 /
  EAS 504 / AFR 661), which bounds the study-heavy split to ~300–400 study rows.

## Risks and Mitigations
- The study-heavy split shrinks the reference panel, which can weaken the AF-LRT
  baseline. Mitigation: report reference $n$ alongside, and keep a partition that
  leaves reference $\ge$ 120.
- The pooled arm may still be PC>AF. Mitigation: publish the measured panel-cap
  bound; the semi-synthetic arm carries the large-$n$ side.
- B3 phased parse time. Mitigation: use the existing `real_phased_chr*` caches.

## Final Checklist
- [x] Main paper tables covered (B1 fills the last real operating point)
- [x] Novelty isolated (no new mechanism; operating-point closure)
- [x] Simplicity defended (no new machinery)
- [x] Frontier contribution justified or explicitly not claimed (no frontier primitive)
- [x] Nice-to-have runs separated from must-run runs (B1 must-run; B2–B4/B6 nice-to-have)
