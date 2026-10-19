# Experiment Plan (Round 13 — a stronger real-data claim)

**Problem**: The paper's strongest real-data statement is currently the *negative*
one — on real 1000G (chr22/1/2, 3 super-populations mixed) the PC channel is
present but *below* AF, so PC-dominance is scoped to the simulated panels. That is
honest but weak: a reviewer reads it as "the headline only lives in the simulator."
The paper's own mechanism claim, however, says the ordering is governed by
*measured structure* (the Fst sweep crosses at measured ≈0.04–0.06). A first
probe on real 1000G confirms the crossover is real and reachable: within a single
super-population (between-population Fst≈0) the study-fitted PC channel reaches
AUC 0.96–1.00, **above** AF (0.81–0.87), gap **+0.13 to +0.15**; across three
super-populations (Fst≈0.086) the ordering flips to AF-first (gap −0.01 to
−0.03). So we can state PC-dominance as a **real-data** phenomenon at the
operating point where it actually applies, instead of a simulated-panel caveat.
**Method Thesis**: The ordering is a function of measured structure, and that
function is now measured *on real genomes* on both sides of its crossover, not
just in the simulator.
**Date**: 2026-10-07

## Claim Map
| Claim | Why It Matters | Minimum Convincing Evidence | Linked Blocks |
|-------|-----------------|-----------------------------|---------------|
| C6 The PC-dominance ordering is a measured-structure law, demonstrated on real genomes on both sides of its crossover | Turns the paper's one negative real-data result into a positive, scoped, real-data law; retires "simulation-only" | Real 1000G, same attack code: PC>AF within a super-population (Fst≈0) and AF>PC across super-populations (Fst≈0.05–0.09), monotone in measured Fst, with the disjoint-negative + structure controls | B1, B2 |
| C7 The within-population PC dominance is not a small-n artefact and holds across independent chromosomes/populations | Blocks "you found it at n=450 on one chromosome" | ≥3 chromosomes × ≥2 within-population cohorts, 5 split seeds each, gap>0 in every cell | B1, B3 |
| C8 The PC channel is attacker-implementable from the release on real within-population data (not just simulated) | The paper's mechanism claim is simulation-scoped; this lifts it | Real within-population study-fitted loadings: full loadings recover, eigenvalue-only / reference-fitted-basis controls at chance | B2 |

## Paper Storyline
- **Main paper must prove**: C6 (the ordering is a real-data, measured-structure
  law). This is the headline real-data upgrade.
- **Appendix can support**: C7 (multi-chromosome/population replication), C8
  (granularity + controls on real within-population data).
- **Experiments intentionally cut**: a real biobank / UKB panel (not obtainable);
  CNV/rare-variant channels (out of scope); a full 22-autosome sweep (chr1/2/22/3
  suffice to show it is not one chromosome).

## Experiment Blocks

### Block 1 (B1): Real-data structure ladder across super-populations
- Claim tested: C6, C7.
- Why it exists: the crossover is currently shown only on the simulator; a
  measured-structure law must be shown where measured structure actually varies —
  real genomes.
- Dataset / split / task: real 1000G phase-3 genotypes, the existing
  `real_multi.get_real_panel_chr`. Cohort definitions from the panel's 26
  `pop` / 5 `super_pop` labels: single-population (GBR, FIN, CHB, JPT, YRI, LWK,
  …), 2-population mixes, and 3/5-super-population mixes. Each cohort: 150 study /
  150 reference / 150 non-member per group, disjoint 3-way split (`three_way_split`),
  5 split seeds.
- Compared systems: AF-LRT (study-fitted AF vs reference AF) vs the study-fitted
  PC projection distance (`common.pca` loadings, study mean), plus the reference-
  fitted-basis and eigenvalue-only controls.
- Metrics: `(af_auc, pc_auc_calibrated, gap=pc-af, fst_like_measured,
  pc_refbasis, pc_eigonly, n_study, M)` per cohort, aggregated over seeds.
- Setup: CPU, project env; chr22 first (fast), then chr1/chr2 from cache.
- Success criterion: PC>AF (gap>0) in every single-population cell and AF>PC in
  every ≥3-super-population cell, with gap monotone-decreasing in measured Fst,
  and both controls at chance in every cell.
- Failure interpretation: if the gap does not flip with measured structure on
  real data, the crossover is simulator-specific and C6 is withdrawn (report the
  real ladder as flat).
- Table / figure target: a new real-data figure (gap vs measured Fst, one point
  per cohort, colored by population count) and a sentence in the abstract.
- Priority: MUST-RUN.

### Block 2 (B2): Real within-population mechanism + granularity
- Claim tested: C6, C8.
- Why it exists: the paper's mechanism (in-sample-ness, loadings-not-eigenvalues)
  is shown on real data only at the *sub-AF* magnitude; we must re-verify it at
  the *dominant* magnitude.
- Dataset / split / task: the within-population cohorts from B1 (GBR and one more),
  chr22 + chr1.
- Compared systems: top-k loadings (k=1,3,5,10) vs eigenvalue-only vs the
  reference-fitted basis vs an out-of-sample (held-out-site) basis.
- Metrics: `(topk_auc[k], eig_only, refbasis, oos_transfer)`.
- Setup: CPU, minutes; reuse `s13_real_loadings`/`s14_real_mechanism` code paths.
- Success criterion: top-k recovers the leak, eig-only and reference-fitted basis
  at chance, out-of-sample transfer at chance — the same pattern as the simulator,
  now at PC>AF magnitude on real data.
- Failure interpretation: if the real within-population leak survives an
  independent basis, the mechanism claim is refuted on real data (a real revision).
- Table / figure target: one mechanism row added to the real-data table.
- Priority: MUST-RUN.

### Block 3 (B3): Independent chromosomes + a second within-population cohort
- Claim tested: C7.
- Why it exists: block B1's within-population result must not be chr22- or
  EUR-specific.
- Dataset / split / task: chr1, chr2 (cached), plus a new chr3 download
  (`ALL.chr3…genotypes.vcf.gz`, ~1.1 GB, one-time parse); within-population
  cohorts GBR, CHB, YRI on each chromosome.
- Compared systems: the B1 cell.
- Metrics: the B1 keys.
- Setup: CPU; chr1/2 cached, chr3 parse ~20 min.
- Success criterion: gap>0 in every within-population × chromosome cell (≥9 cells),
  5 seeds each, and the controls at chance.
- Failure interpretation: a chromosome where the gap vanishes bounds the claim to
  the chromosomes where it holds (reported per-chromosome).
- Table / figure target: appendix replication table.
- Priority: MUST-RUN.

## Run Order and Milestones
| Milestone | Goal | Runs | Decision Gate | Cost | Risk |
|-----------|------|------|---------------|------|------|
| M0 | sanity: one single-pop and one 3-super-pop cell reproduce the probe | B1 2 cells | gap>0 vs gap<0 confirmed with controls at chance | <10 min CPU | split differs from probe |
| M1 | real structure ladder (chr22) | B1 full | monotone gap in measured Fst | ~30 min CPU | population-pair scarcity |
| M2 | real within-pop mechanism | B2 | loadings recover, controls at chance | ~15 min CPU | real OOS differs from sim |
| M3 | multi-chromosome within-pop | B3 (chr3 parse + sweep) | gap>0 in all ≥9 cells | ~1.5 h CPU | chr3 download |
| M4 | run the paper-review loop | — | reviewer accepts the real-data law | — | reviewer asks for more pops |

## Compute and Data Budget
- Total: CPU-only; ~2.5 h wall. No GPU needed (panels are 1.6k individuals).
- Data: chr3 1000G VCF (~1.1 GB) one-time download + parse; chr1/2/22 cached.
- Human evaluation: none.
- Biggest bottleneck: the per-cell `np.linalg.svd` on ~450×50k — already done in
  `s7_real`/`s16`, minutes per cell.

## Risks and Mitigations
- A single-population cohort may have too few samples after the 3-way split: use
  26 `pop` labels (each 60–113 samples) so per-group n is smaller but the split is
  valid; report n per cell.
- The within-population "no between-population structure" claim must be measured,
  not asserted: report measured `fst_like` (≈0 within, ≈0.05–0.09 across) per cell,
  and the reference-fitted basis as the structure control.
- Reviewer may say the within-population operating point is contrived: frame it as
  the *realistic* fine-mapping / population-specific-study case, with the across-
  population case as the boundary, both measured.
- chr3 download fails: fall back to chr1+chr2 (already 2 independent chromosomes).

## Final Checklist
- [ ] Real-data gap flips with measured Fst, controls at chance (C6)
- [ ] Within-population dominance replicated across ≥3 chromosomes / ≥2 pops (C7)
- [ ] Real within-pop mechanism + granularity reproduced (C8)
- [ ] The across-population (AF-first) boundary retained
- [ ] Nice-to-have separated from must-run
