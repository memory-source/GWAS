# Experiment Plan — Fair Comparison and Interface Confirmation (round 22, 2026-10-17)

**Center claim under test** (narrowed from the round-19 formulation): a cooperating
analysis partner that keeps every per-SNP genotype count unchanged, and makes only
bounded input modifications, obtains membership information beyond the strongest
passive adversary through a small number of bounded-precision final PCA outputs; the
capability has interpretable precision, scale, and permission boundaries.

**Why this round exists.** The round-19/20/21 results established the identity, the
marginal-preservation lemma, the coverage growth, and an accumulated-output increment
measured against the *larger of two separate* passive baselines. The round-22 GPT
analysis raises two questions that decide whether the route continues:

1. **Fair comparison.** The margin was measured against `max(af_cond, spec_cond)`, the
   better of two *single* channels. A real adversary holds both channels and fuses
   them. If the increment survives only against the stronger single channel and
   vanishes against the fused channel, the claim must shrink.
2. **Interface confirmation.** The attack assumes the attacker can obtain the joint
   top-k loadings and eigenvalues, the public center/scale, and the right to re-run on
   the same honest cohort after modifying its own input, at the printed precision. This
   must be established item by item, and the attack must survive an *approximate*
   (randomized) solver, not only the exact eigen-decomposition.

**Staging rule (from the GPT analysis).** Run the fair-comparison and interface blocks
first. Only if they pass do the precision-robustness, known-member-fraction, independent
-cohort and science-cost blocks proceed at scale. This plan lists all blocks; the first
two gate the rest.

## Claim Map

| Claim | Why it matters | Minimum convincing evidence | Blocks |
|-------|----------------|-----------------------------|--------|
| J1 — the accumulated-output increment survives the *fused* passive adversary (conditioned AF joined with the conditioned single release, weights fitted on the calibration half) | decides whether the active mechanism adds anything over a strong passive adversary | paired `gap_vs_joint_passive` with 95% CI at q=4/q=8, on both chromosomes, per split; reported in AUC *and* TPR@5%FPR | B67 |
| J2 — the capability is not preprocessing information: the public-mean/scale-only channel and the "same input, repeated outputs" control do not carry it | rules out that the leak is the center/scale or repeated sampling | gap of the q-attack over (i) the center/scale-only score and (ii) the passive repeats | B67 |
| J3 — the interface is reproducible: per item, the attacker's access to loadings, eigenvalues, center/scale, re-run right, and printed precision is stated and demonstrated on a real software export | the threat model must be real, not assumed | machine-readable interface manifest + a PLINK 2.0 `--pca allele-wts` read-back (exact) | B68 |
| J4 — the attack survives an approximate (randomized) solver | real deployments use `--pca approx`; an exact-equation attack may not transfer | identity residual + stacked-attack gap under `--pca approx`; span growth from solver seed changes on the SAME input | B68 |
| J5 — the reliable-direction set is bounded by the export precision, not only by the fixed 1e-9 rank cutoff | a direction that is numerically new is not necessarily usable under rounding | per-direction reliability floor derived from the printed precision; coverage/membership recomputed on the reliable set at float64/float32/6dp/4dp | B69 |
| J6 — the increment is characterized against the known-member fraction (0/5/10/50%) | the calibrated-attacker advantage must be stated for the cooperation scenarios it holds in | `gap_vs_joint_passive` at each known-member fraction; 0% oriented on shadow data | B70 |
| J7 — under protection the honest GWAS keeps usable power at a reported cost | a defense is only interesting if the science survives; lambda-GC and min-p alone are insufficient | release-strategy table: FP rate, power, effect bias, individual QC pass, per-release cost | B71 |
| J8 — independent-subject replication | another chromosome or a re-split is not an independent cohort | HGDP availability assessment + duplicate/relatedness/ancestry checks, then a frozen-method replication | B72 |

## Experiment Blocks

### Block 67 (B67): the fair joint passive adversary — MUST-RUN, GATE
- Claim J1, J2. Script `code/s67_joint_passive.py` →
  `results/raw/s67_joint_passive_chr<C>.json`.
- For each cohort (`within:EUR`, `super:3pop`, `within:EAS`) × 5 split seeds:
  - build the two conditioned passive channels on the SAME held-out split as the
    spectral attack: `af_cond` (composition-balanced removal of the known calibration
    members from the released AF) and `spec_cond` (the conditioned single release);
  - **fuse** them: z-score both channels on the *calibration half*, fit a logistic
    regression of the known member labels on the fused features (the attacker knows
    the calibration members), and evaluate the fused score on the held-out half;
  - the joint adversary is the better of the fused logistic score and the stronger
    single channel on the calibration half (still oriented/evaluated held-out);
  - run the q=1/2/4/8 random and adaptive stacks exactly as in `s61.cell`, then report
    `gap_vs_joint_passive` in AUC and TPR@5%FPR with bootstrap CIs and the positive
    split count.
- Controls in the same cell: the center/scale-only score (no spectral output) and the
  passive-repeat score (same input, q outputs).
- Decision rule: if the q=8 `gap_vs_joint_passive` CI excludes zero on both
  chromosomes in at least the cross-population cohort, J1 holds and the route
  continues; otherwise the claim is narrowed to the fused-baseline-wins statement.

### Block 68 (B68): interface manifest + approximate PCA — MUST-RUN, GATE
- Claim J3, J4. Script `code/s68_interface.py` →
  `results/raw/s68_interface_chr<C>.json`.
- Interface manifest: a JSON object listing, per item, the attacker's access and the
  artifact that demonstrates it: joint top-k loadings (`--pca allele-wts` →
  `.eigenvec.allele`), eigenvalues (`.eigenval`), public center/scale (the released
  pooled AF), the right to re-run on the same honest cohort (the joint VCF is
  rewritten per release), and the printed precision (6 significant digits). Each entry
  carries `{obtainable: bool, demonstrated_by: path, note}`.
- Approximate solver: re-run the joint PCA with `--pca <k> allele-wts approx` (the
  randomized Galinsky et al. algorithm) on the same joint VCF; read back; compute the
  identity residual `max‖C_H V_t − Y_t‖ / max‖C_H V_t‖` and the q8 stacked attack vs
  the conditioned single release. Compare with the exact-solver residual.
- Solver randomness on a FIXED input: run `--pca approx` with two different `--seed`
  values on the SAME joint VCF; measure how far the returned top-k span moves and
  whether stacking the two spans (a "same-input, repeated output" attack) grows the
  span. This is the solver-randomness control for J2/J4.

### Block 69 (B69): precision-reliable directions — CONDITIONAL on the gate
- Claim J5. Script `code/s69_reliable_dirs.py` →
  `results/raw/s69_reliable_dirs_chr<C>.json`.
- Derive a per-direction reliability floor from the printed precision: for each
  precision (float64, float32, 6 dp, 4 dp) estimate the operator perturbation norm
  `‖Δ(C_H V − Y)‖` and mark direction `i` usable only if its singular value in the
  stacked `Q` (or the released `V_t`) exceeds the perturbation induced at that
  precision. Recompute visible share, scatter reconstruction, and membership on the
  reliable set only, and compare with the fixed 1e-9 cutoff.

### Block 70 (B70): known-member-fraction sweep — CONDITIONAL on the gate
- Claim J6. Script `code/s70_known_fraction.py` →
  `results/raw/s70_known_fraction_chr<C>.json`.
- Sweep the fraction of honest members whose labels the attacker knows
  {0%, 5%, 10%, 50%}. At 0% the orientation and the conditioned baselines are fixed on
  the independent shadow data. Report `gap_vs_joint_passive` at each fraction.

### Block 71 (B71): scientific cost under protection — CONDITIONAL on the gate
- Claim J7. Script `code/s71_protection_cost.py` →
  `results/raw/s71_protection_cost_chr<C>.json`.
- Compare release strategies (a fixed single pass; an input-varying pass; an
  input-consistency-constrained pass) on the honest analysis: GWAS false-positive rate
  at a fixed threshold, power at the planted causal markers, causal effect-size bias,
  per-individual QC pass rate (heterozygosity / missingness), and per-release cost —
  beyond lambda-GC and the minimum p-value.

### Block 72 (B72): independent-subject replication — CONDITIONAL, last
- Claim J8. Assess HGDP availability (Sanger/CEPH public panel), run duplicate /
  relatedness / ancestry-matching checks, then replicate the frozen method on an
  independent cohort; use same-source members and non-members only. Fixed threshold on
  the calibration set; report test-set correct / false-positive counts with intervals.

## Run Order
B67 → B68 (gates). If both pass: B69 → B70 → B71 → B72. Reuse `code/common.py`,
`code/s61_spectral_tomography.py` (`split_roles`, `zscore`, `joint_baseline`, `topk`,
`swap_block`, `_membership_auc_from_operator`, `af_baselines`, `cell`),
`code/s48_real_structure.py` (`load_for`, `_orient_auroc`),
`code/s51_real_lowstruct.py` (`_sup_label`), `code/s63_plink_export.py` (VCF write,
PLINK run, read-back).

## Compute Budget
CPU only. Real 1000G chr22/chr1 prefixes (1.5 kb thinning); honest cells ≤ 450 rows;
one dense eigen-decomposition per release; a PLINK 2.0 run per release for B68.
Project env `/data/lvvyongji1/GWAS/env/bin/python`; PLINK static binary
`scratch/plink2bin/plink2`.

## Scope Words
"real genotype driver + mechanism-clear simulated phenotype"; no "comprehensive",
"robust", or "extensive". Every range enumerated from the JSON. The increment is
stated against the fused passive adversary, not the stronger single channel.
