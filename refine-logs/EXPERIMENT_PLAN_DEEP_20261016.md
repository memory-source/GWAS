# Experiment Plan — Selection-Input Spectral Tomography (deep-search round)

**Problem**: The passive paper measures leakage from a *fixed* release. The
deep-search direction (`research-plan/20261016_deep_search/DEEP_DIRECTIONS_ZH.md`)
proposes the active question: when a cooperating analysis partner may vary its
own AF-preserving input across re-fits of a joint PCA and reads only the final
truncated output each time, how many independent measurements of the honest
cohort can it design? The synthetic algebra witness (`algebra_witnesses.py`)
shows the mechanism exists in an artificial independent-locus model. Nothing on
real genomes, no control baseline, no adaptive selection has been run.

**Method Thesis**: Model the joint PCA final output as a *measurement oracle*:
with honest scatter $C_H$ and the partner's known scatter $K_t$, the released
top-$k$ eigenvectors $V_t$ and eigenvalues $\Lambda_t$ give
$Y_t = V_t\Lambda_t - K_t V_t = C_H V_t$. Stacking $q$ outputs,
$C_H Q = Y$ with $Q=[V_1,\dots,V_q]$, so $\hat C_H = Y Q^{+}$ whenever $Q$ has
full row rank. The scientific question is how large a fraction of $C_H$ a small
number of *real, QC-passing, AF-preserving* block swaps can cover, and whether
that coverage turns into a membership/family advantage over the strongest
passive AF baseline.
**Date**: 2026-10-16
**Design constraint from the direction doc**: real genotype + mechanism-clear
*simulated* phenotype; adaptive selection may read only partner data, past
outputs, and independent shadow data — never the honest private matrix.

## Claim Map
| Claim | Why It Matters | Minimum Convincing Evidence | Linked Blocks |
|-------|-----------------|-----------------------------|---------------|
| (A) The final-output identity $Y_t=C_H V_t$ holds exactly under the released preprocessing | Without it the "measurement" reading is wrong | Max residual over all releases/cells below ~1e-9 (float64) | B0 |
| (B) Each AF-preserving block swap leaves per-SNP genotype histograms, AF, HWE and per-SNP missing counts unchanged | This is the threat-model premise (data passes marginal QC) | Exact histogram equality per release + AF/variance max-abs-diff | B0 |
| (C) Visible-subspace coverage and target-scatter recovery grow with the number of final outputs, and adaptive block choice reaches a target orientation faster than random | The core new capability | subspace-share + relative $\hat C_H$ error + target-direction residual vs #outputs ∈ {1,2,4,8}; adaptive vs random at matched outputs | B1, B2 |
| (D) The accumulated measurement yields a membership/family-identification increment over the strongest passive AF baseline | Turns a mechanism into a privacy capability | held-out-half AUROC of an $\hat C_H$-based projection-distance vs AF-LRT, matched information, disjoint negatives, CI + permutation null | B3 |
| (E) The per-reanalysis scientific cost is small at the useful coverage levels | The route is only interesting if the honest analysis survives | per-release $\lambda_{GC}$, top-PC score correlation, GWAS power on a planted signal | B4 |
| (F) AF-preservation is not a full-authenticity guarantee | Scope honesty | same-partner-real-data durability probe: within-SNP histogram kept, but individual heterozygosity/relatedness statistics shift | B5 |

## Paper Storyline
- **Main paper** (security/ML venue): the identity (A), the marginal-preservation
  lemma (B), the coverage/recovery curves (C), and the privacy increment over a
  strong AF baseline (D). One headline figure = the "most decisive next figure"
  of the direction doc.
- **Appendix**: the synthetic algebraic witness re-run under the same interface
  (B0) plus the science-cost detail (E) and durability probe (F).
- **Cut**: the blind rare-mixture route (high-risk, separate paper) and the
  transferable-structure route (statistical-genetics paper). Not merged here.

## Experiment Blocks

### Block 0 (B0): interface witness on real geometry + preservation audit — MUST-RUN
- Claims (A),(B). Meaningful on the real marker layout so the later real cells are
  interpretable, beyond the 24-SNP toy.
- Build: honest scatter from real 1000G chr22 with real LD; partner matrix =
  a disjoint 1000G draw (so LD is real). Input change = swap one LD-consistent
  SNP block between two partner rows.
- Checks: (i) max $|(C_H V_t) - Y_t|$ over releases and seeds; (ii) per-SNP
  genotype histogram identical to the initial partner matrix; (iii) AF, variance,
  HWE and per-SNP missing counts identical.
- Script: `code/s61_spectral_tomography.py` (`witness_realcell`).
- Output: `results/raw/s61_spectral_tomography.json` (`witness` + `preservation`).
- Priority: MUST-RUN.

### Block 1 (B1): subspace coverage vs number of outputs — MUST-RUN
- Claim (C). Fixed honest cohort, x-axis #final outputs {1,2,4,8}.
- Measure: rank of $Q$; share of $C_H$'s squared spectrum captured by the
  measurement span; relative Frobenius error of $\hat C_H=Y Q^{+}$; smallest
  reported singular value of $Q$ (conditioning).
- Strategies: passive (same input each release) / random AF-preserving block
  swap / adaptive (next block chosen to reduce the residual of a named target
  direction, using only partner data + past outputs + independent shadow).
- Script: main runner over strategies × #outputs × 5 seeds.
- Priority: MUST-RUN.

### Block 2 (B2): target-direction recovery — MUST-RUN
- Claim (C). For a named candidate direction $x$ (e.g. an ancestry direction read
  from independent shadow data), report residual $\|x - Qa\|$ and the
  $|x^\top C_H x - x^\top Y a|$ error vs #outputs, adaptive vs random.
- Priority: MUST-RUN.

### Block 3 (B3): membership / family increment over the strong AF baseline — MUST-RUN
- Claim (D). Disjoint 3-way split (honest members / attacker reference / non-member
  targets). Score targets by projection distance in the $\hat C_H$-recovered
  subspace (and/or by a member-direction likelihood), held-out-half oriented,
  against: unconditioned AF-LRT, conditioned AF (balanced removal of known
  calibration members), and passive single-release PC.
- Report per-split gaps + bootstrap CI + stratified permutation null (reuse
  `s48/_orient_auroc`, `s54` machinery).
- Priority: MUST-RUN.

### Block 4 (B4): scientific-cost side — NICE-TO-HAVE
- Claim (E). Per release: $\lambda_{GC}$ (corrected), correlation of the honest
  top-PC scores across re-fits, and GWAS power on a planted simulated causal
  signal via `common.marginal_association`.
- Priority: NICE-TO-HAVE (report; does not gate).

### Block 5 (B5): durability / scope probe — NICE-TO-HAVE
- Claim (F). Same swap sequence applied to a *real* partner sub-panel (no
  synthetic AF guarantee) to show the within-SNP histogram check is what is
  preserved while individual-level structure statistics move.
- Priority: NICE-TO-HAVE.

## Run Order
B0 → B1 → B2 → B3 (all when the runner is up), then B4/B5.
Reuse: `code/common.py` (AF-LRT, pca, marginal_association, auroc,
stratified_permutation_auroc), `code/real_multi.py` (1000G chr22 loader),
`code/s48_real_structure.py` (`_orient_auroc`, real cohort loading),
`code/s54_real_permnull.py` (permutation null + CI helpers).

## Compute Budget
CPU only; real chr22 panel ≤ 20k markers thinned to ~1.5 kb spacing; honest cells
n ≤ ~450. Estimator cost is a small dense eigen-decomposition per release.
Runs on the project env `/data/lvvyongji1/GWAS/env/bin/python`.

## Scope Words
"real genotype driver + simulated phenotype" everywhere; no "comprehensive",
"robust", or "extensive". Every reported range enumerated from the JSON.
