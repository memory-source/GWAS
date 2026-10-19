# Experiment Results — Fair Comparison and Interface Confirmation (round 22, 2026-10-17)

Direction: `research-plan/20261016_deep_search/DEEP_DIRECTIONS_ZH.md`.
Plan: `refine-logs/EXPERIMENT_PLAN_ARB_20261017.md`.
New scripts: `code/s67_joint_passive.py`, `code/s68_interface.py`,
`code/s69_reliable_dirs.py`, `code/s70_known_fraction.py`,
`code/s71_protection_cost.py`, `code/s72_printed_matched.py`.
Raw: `results/raw/s67_joint_passive_chr{22,1}.json`, `s68_interface_chr22.json`,
`s69_reliable_dirs_chr{22,1}.json`, `s70_known_fraction_chr{22,1}.json`,
`s71_protection_cost_chr{22,1}.json`, `s72_printed_matched_chr{22,1}.json`.
Reused: `code/s61_spectral_tomography.py` (split roles, z-score, baseline, top-k,
block swap, membership), `code/s63_plink_export.py` (VCF write, PLINK run, read-back),
`code/s48_real_structure.py`, `code/s51_real_lowstruct.py`, `code/common.py`.

All cells: real 1000G chr22/chr1 prefixes (1.5 kb thinning), 5 split seeds per cohort,
held-out-half orientation, the same 5-role disjoint split as round 19.

## J1 — the increment survives the FUSED passive adversary (the decisive comparison)

The round-19/20/21 margin was measured against the stronger of two SINGLE passive
channels. Here the passive adversary is given the strongest of FOUR candidates —
`af_cond`, `spec_cond` (the two single conditioned channels), `fused_z` (their
label-free z-sum, orientation = mechanism sign), and `fused_cal` (a logistic fusion
whose weights are fitted on the unconditioned calibration half) — and the bar is the
best of them on the HELD-OUT half (`oracle`), i.e. the strongest passive adversary
these constructions can build. Using an oracle bar is conservative against the active
claim: no candidate-selection rule could be stronger. The AUC bar is the AUROC-best
candidate and the TPR bar is the TPR-best candidate (the two are reported separately,
since the strongest AUC passive need not be the strongest TPR passive). The active
attacker is the q = 1/2/4/8 random/adaptive stack, scored from the recovered operator;
`gapA` = the paired q8 AUC increment over the AUC bar, `gapT` = the paired q8 TPR@5%FPR
increment over the TPR bar, 10k-bootstrap 95% CI, 5 splits.

### chr22 (`results/raw/s67_joint_passive_chr22.json`)
| cohort | af_cond | spec_cond | fused_z | fused_cal | AUC bar | random q8 gapA | TPR bar | q8 gapT | q8 posA/posT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| within:EUR | 0.882 | 0.985 | 0.989 | 0.993 | 0.993 | +0.0066 [+0.0027,+0.0111] | 0.960 | +0.040 [+0.017,+0.066] | 5/5 |
| super:3pop | 0.784 | 0.757 | 0.840 | 0.837 | 0.842 | +0.0656 [+0.0541,+0.0822] | 0.436 | +0.144 [+0.077,+0.210] | 5/5 |
| within:EAS | 0.915 | 0.987 | 0.994 | 0.996 | 0.996 | +0.0037 [+0.0011,+0.0061] | 0.983 | +0.017 [+0.007,+0.026] | 4/5 |

### chr1 (`results/raw/s67_joint_passive_chr1.json`)
| cohort | af_cond | spec_cond | fused_z | fused_cal | AUC bar | random q8 gapA | TPR bar | q8 gapT | q8 posA/posT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| within:EUR | 0.906 | 0.965 | 0.984 | 0.985 | 0.987 | +0.0130 [+0.0093,+0.0181] | 0.926 | +0.074 [+0.056,+0.096] | 5/5 |
| super:3pop | 0.805 | 0.748 | 0.846 | 0.839 | 0.850 | +0.0443 [+0.0272,+0.0599] | 0.460 | +0.015 [−0.063,+0.092] | 5/3 |
| within:EAS | 0.933 | 0.980 | 0.991 | 0.991 | 0.992 | +0.0048 [+0.0022,+0.0073] | 0.980 | +0.017 [+0.003,+0.030] | 5/3 |

(`gapA` = AUC increment over the AUC bar; `gapT` = TPR@5%FPR increment over the
metric-specific TPR bar; `posA/posT` = positive AUC-gap / TPR-gap splits. The AUC bar is
`fused_cal` on 4/5 chr22 EUR splits and `fused_z`/`fused_cal` mixed elsewhere; the TPR
bar is `fused_z`/`fused_cal` mixed. Every CI in this section is a 10k percentile
bootstrap over the five repeated splits of ONE fixed panel — conditional split
uncertainty, not independent-subject generalization.)

The accumulated q8 attack beats the metric-specific passive bar in all three cohorts on
both chromosomes with a positive AUC CI; 5/5 AUC-positive splits except chr22 EAS (4/5).
At the practical TPR@5%FPR reading the increment is measured against the strongest TPR
passive bar and is larger (EUR +0.040/+0.074, 3pop +0.144/+0.015, EAS +0.017/+0.017).
The TPR gaps are positive in 5/5 splits for chr22 EUR/3pop and chr1 EUR, 4/5 for chr22
EAS, and 3/5 for chr1 EAS/3pop; the conditional bootstrap interval for the mean excludes
zero everywhere except chr1 3pop (chr1 EAS mean interval [+0.0033,+0.0298]). At q4 the positive
AUC-gap CI excludes zero only for chr1 within:EUR (adaptive +0.009 [+0.005,+0.014] 5/5);
the chr22 EAS random and both chr1 3pop q4 cells have a CI entirely below zero (chr22 EAS
adaptive q4 crosses zero). The q8 checkpoint is where the effect is reliable in every
cohort. The passive strategy never beats the bar at any q.

## J2 — the center/scale-only channel and the same-input repeat do not carry it

- Center/scale-only (`af_pooled`, the public mean/scale with no per-member release):
  chr22 0.791 / 0.695 / 0.812, chr1 0.783 / 0.688 / 0.810 — well below the joint
  passive adversary and the active stack. The effect is not the public preprocessing.
- Same-input repeat: the passive strategy's q1..q8 rows are identical by construction
  (one fixed input), so repeated releases of the same input add no measurement.
- Same-input solver randomness (s68): PLINK `--pca approx` with two different `--seed`
  values on the SAME joint VCF moves the top-k span by a mean per-split maximum-angle
  sine of 6.0e-8 (across-split max 8.2e-8; reported as the mean); stacking the two
  seed spans gives a TPR@5%FPR gain of 0.000. Repeated solver runs on a fixed input
  are not a new measurement.

## J3 — the interface is reproducible (machine-readable manifest)

`s68_interface_chr22.json` `interface_manifest` states each item. The manifest
demonstrates LOCAL SOFTWARE CAPABILITY under an ASSUMED access model: the loadings and
eigenvalues items are backed by real PLINK artifacts, and the re-run/center-scale items
are what the attack constructs locally from the plan's granted information. Whether a
deployed federated pipeline authorizes a partner to re-run against an unchanged private
cohort, and whether it publishes an AF artifact, are threat-model assumptions, not
findings. The honest private matrix is never read:

| item | obtainable | demonstrated by |
|---|---|---|
| joint top-k loadings | yes | PLINK 2.0 `--pca <k> allele-wts` → `.eigenvec.allele` (ALT rows) |
| eigenvalues | yes | `.eigenval` (= scatter-eigenvalues / marker count) |
| public center/scale | yes | the released pooled AF (honest+partner mean) |
| re-run right | yes | the joint VCF is rewritten with the partner's modified input per release |
| printed precision | yes | the `.eigenvec.allele` decimals (6 significant digits) |
| honest private matrix | no | the attack never reads `C_H` or the honest genotypes |

## J4 — the attack survives the approximate solver

PLINK `--pca approx` (the randomized algorithm; its log records Krylov-matrix SVD and
range recovery, so it is genuinely executed) on the same joint VCFs gives the SAME
identity residual (6.242e-6 max in both) and the SAME q8 stacked AUROC (0.998637 in
both). An independent export comparison found the eigenvalues identical (difference
0.0) and the spans near-identical; the eigenvector columns differ only by sign flips on
near-degenerate eigenpairs (e.g. PC8 entries +3.046 vs −3.046), so after aligning
signs the maximum entry difference across all 40 pairs is about 1e-6 — a sign/rotation
artifact, not a span change. The matched printed-output comparison (J9) reruns the
whole attack through `approx` and reproduces the same conclusion. Caveat: PLINK's
`approx` targets n > 5k; at n = 240 here it introduces no detectable span change, so
this is a lower bound on the approx stress, not a biobank-scale proof.

## J5 — heuristic retained rank under simulated quantization (offline float64 reference)

`s69_reliable_dirs`: the attack (s61.metrics) scores membership by projection onto the
range of the recovered honest-cohort operand `Y_t = V_t Λ_t − K_t V_t = C_H V_t`
(stacked), which is what the released V/Λ actually expose. With the fixed 1e-9 cutoff
the reported rank RISES under rounding (it counts spurious near-zero directions).
Applying a reliability floor derived from the measured export perturbation
`dY = ‖Y_prec − Y_f64‖ / ‖Y_f64‖` removes them; the reliable set is re-measured for
coverage and membership. This is a descriptive count under one plausible floor, not a
certified usable-rank theorem: the perturbation is measured against the float64 operand,
which the attacker never sees, so the floor is an offline construct for the measurement
rather than a rule the deployed attack could compute from the printed precision alone.

| cohort | precision | geometric rank | reliable rank | dY | visible_share (reliable) | AUROC (geo → reliable) |
|---|---|---:|---:|---:|---:|---:|
| within:EUR | float64 | 65 | 65 | 0 | 0.941 | 1.000 → 1.000 |
| within:EUR | float32 | 78 | 53 | 6.7e-8 | 0.913 | 1.000 → 1.000 |
| within:EUR | 6 dp | 78 | 30 | 3.2e-5 | 0.830 | 0.996 → 0.993 |
| within:EUR | 4 dp | 78 | 11 | 3.2e-3 | 0.677 | 0.963 → 0.920 |
| super:3pop | float64 | 61 | 61 | 0 | 0.968 | 0.907 → 0.907 |
| super:3pop | 6 dp | 80 | 16 | 3.3e-5 | 0.938 | 0.806 → 0.711 |
| super:3pop | 4 dp | 80 | 10 | 3.3e-3 | 0.929 | 0.678 → 0.655 |
| within:EAS | 6 dp | 78 | 29 | 2.8e-5 | 0.834 | 0.990 → 0.990 |
| within:EAS | 4 dp | 78 | 11 | 2.8e-3 | 0.686 | 0.953 → 0.923 |

(chr1 reproduces the same ordering: EUR 6 dp 74→27, 4 dp 74→11; 3pop 6 dp 80→17,
4 dp 80→10.) The quantization is decimal-place rounding of the unit-norm eigenvector
entries `np.round(V, 6)` / `np.round(V, 4)`, a SIMULATED quantization, not PLINK's
six-significant-digit print of the differently scaled allele weights — the two are not
interchangeable, so these rows bound the sensitivity of the attack to eigenvector
rounding rather than certify a PLINK-specific precision ceiling. The real printed-export
precision is measured separately by s63 (exact-vs-approx read-back of actual PLINK
outputs) and by the matched printed cell (J9). The geometric cutoff is precision-blind;
the reliable set shrinks as the precision falls. At 6 dp the within-population reliable
AUROC barely moves (EUR 0.996→0.993, EAS 0.990→0.990); at 4 dp the within-population
reliable AUROC falls visibly (EUR 0.963→0.920, EAS 0.953→0.923). The cross-population
3pop reliable AUROC drops more sharply still (0.806→0.711 at 6 dp), so there the removed
directions carried membership signal. The rank numbers are outputs of this
perturbation-based reliability heuristic on the operand stack, not a certified
usable-rank theorem.

## J6 — the increment holds down to 0% known members

`s70_known_fraction`: the fraction is a share of the HONEST COHORT `N` (the drawn
members are restricted to the calibration half, the labelled part), so at nominal 50%
the attacker knows 56–60 of the 120 EUR/EAS members and 178–180 of the 360 3pop
members. The known members are removed (composition-balanced, equal quota per group)
from the released AF and the single release; orientation is the MECHANISM sign (+1) for
both channels at EVERY fraction (a member's genotype matches the study AF better, and
its scatter made it project further along the released range), so no member label is
ever used to fix a sign. The passive bar is the `fused_z` label-free fusion (reachable
at every fraction, no label-fitted weights). The increment shrinks as the passive bar
strengthens with more removable members but stays positive (5/5, CI excluding zero) at
every fraction including 0%:

| cohort | 0% joint / gap | 5% joint / gap | 10% joint / gap | 50% joint / gap |
|---|---:|---:|---:|---:|
| within:EUR | 0.936 / +0.064 | 0.941 / +0.058 | 0.946 / +0.054 | 0.989 / +0.011 |
| super:3pop | 0.724 / +0.184 | 0.729 / +0.178 | 0.740 / +0.167 | 0.842 / +0.065 |
| within:EAS | 0.940 / +0.059 | 0.943 / +0.057 | 0.953 / +0.046 | 0.994 / +0.006 |

(chr1 reproduces; e.g. EUR 0% +0.076 to 50% +0.017, 3pop 0% +0.175 to 50% +0.052.)
The calibrated advantage is therefore NOT load-bearing: the increment is largest when
the attacker knows nothing and shrinks only when the passive adversary is handed the
member labels (the conditioned `fused_z` bar rises to 0.989/0.994 within a population
at 50%), so the within-population headline margin rests at 50% on a near-ceiling gap
while the cross-population 3pop arm keeps a larger margin (+0.065).

## J7 — scientific cost under protection (scoped)

`s71_protection_cost`: honest GWAS (simulated phenotype, 8 causal markers with exposed
per-marker effects `~N(0,0.35)`, fixed p < 1e-4 threshold, ONE QC marker universe for
every arm) under a fixed single pass and under an input-varying pass. The effect
comparison is `effect_agreement` = Pearson r(estimated effect, true effect) over the
causal markers (scale-free), chr22 and chr1:

| cohort | arm | lambda_GC | power (causal) | recovered (mean) | off-causal fraction | effect_agreement |
|---|---|---:|---:|---:|---:|---:|
| within:EUR | honest alone | 1.072 | 0.025 | 0.2/8 | — | — |
| within:EUR | passive / perturbed | 1.095 / 1.096 | 0.025 / 0.025 | 0.2/8 | 5.5e-5 | 0.584 / 0.584 |
| super:3pop | honest alone | 1.136 | 0.050 | 0.4/8 | — | — |
| super:3pop | passive / perturbed | 1.149 / 1.148 | 0.050 / 0.050 | 0.4/8 | 1.9e-4 | 0.833 / 0.833 |
| within:EUR (chr1) | passive / perturbed | 1.158 / 1.158 | 0.000 / 0.000 | 0/8 | 0 | 0.377 / 0.377 |
| super:3pop (chr1) | passive / perturbed | 1.072 / 1.072 | 0.025 / 0.025 | 0.2/8 | 2.6e-4 | 0.872 / 0.872 |

The honest analysis is unchanged between the fixed pass and the perturbed pass (lambda
within 0.001, power and recovered-causal count identical, off-causal discovery
fraction identical, effect_agreement within 0.001). The mean recovered count is
expressed per 8 causal markers (e.g. 0.2/8 = 0.2 recovered out of 8 on average).
`off-causal fraction` is the share of tested markers that were NOT planted as causal and
fell below the 1e-4 threshold. It is a discovery count, not a clean null
false-positive rate: under linkage disequilibrium some off-causal markers are proxies of
causal ones, so a discovery there is not necessarily an error. The chr22 3pop value is
roughly 0.6 unplanted markers per split of ~3100 tested. The marginal partner quantities
are exactly preserved (per-SNP genotype histograms, AF, variance, and the per-marker
scatter), so the released per-SNP counts are unchanged; individual-level statistics do
change slightly (row heterozygosity moves by up to 4.5e-3 across the release sequence,
and the maximum per-SNP variance change is at machine precision, 2.2e-16 to 6e-16), which
is the physical signature of reorderings inside a fixed per-SNP histogram. The honest rows are untouched.
Power is low because the honest
study is small (n = 120/360) with a simulated phenotype; this is a summary-level
stability statement at one study size, not a power curve. The planned
input-consistency-constrained arm and a per-release wall-clock cost were not
implemented.
## J9 — the matched printed-output comparison (the bounded-precision conjunction)

`s72_printed_matched` closes the round-22 review's demand: the earlier J1/J4/J5 cells
measured the accumulated attack (J1) in float64, showed the approx solver leaves the
span unchanged (J4) and bounded the reported direction count (J5), but no single cell
measured the accumulated attack against the passive baseline in the SAME PRINTED OUTPUT.
Here each cohort x split runs q = 8 joint re-fits through real PLINK 2.0
(`--pca 10 allele-wts approx`, `--seed 1000+t`), reads the printed allele weights and
eigenvalues back, forms the attack operand `Y_t = V_t Λ_t − K_t V_t` from the READ-BACK
V/Λ (PLINK's scatter/M normalization) with the partner's own scatter `K_t` known
exactly, and compares the q8 stack against the passive constructions built from the
same printed outputs. The passive comparator is built under the SAME information model as
s67's `spec_cond`: it sees the same released V/Λ and the attack operand `Y_t` (the
conditioned single release subtracts the known calibration members' projection `Z_k V_t`),
fits the calibrated-fusion weights on the UNCONDITIONED calibration features and applies
them to the conditioned channels, and fuses the per-release scalar channels over q (`oracle_q`
is a late-fusion diagnostic over the four specified constructions, not a general statement
about arbitrary processing of the release transcript). Two matched bars are reported —
`oracle_single` (best of the four constructions of the t=0 output alone) and `oracle_q`
(best of the four constructions fused over the q printed outputs).

| cohort | chr | oracle_single | oracle_q | q8 | gap vs oracle_single | gap vs oracle_q | pos (single / q) |
|---|---:|---:|---:|---:|---:|---:|---:|
| within:EUR | 22 | 0.993 | 0.993 | 0.999 | +0.0054 [−0.0003,+0.0110] | +0.0055 [−0.0001,+0.0110] | 4/5 / 4/5 |
| super:3pop | 22 | 0.842 | 0.842 | 0.850 | +0.0082 [−0.0084,+0.0287] | +0.0082 [−0.0084,+0.0287] | 4/5 / 4/5 |
| within:EUR | 1 | 0.987 | 0.987 | 0.997 | +0.0100 [+0.0038,+0.0165] | +0.0100 [+0.0037,+0.0162] | 4/5 / 4/5 |
| super:3pop | 1 | 0.850 | 0.850 | 0.856 | +0.0060 [−0.0149,+0.0237] | +0.0060 [−0.0149,+0.0236] | 3/5 / 3/5 |

The matched printed-output increment is positive in every cell (EUR chr22 +0.0054, 3pop
chr22 +0.0082, EUR chr1 +0.0100, 3pop chr1 +0.0060 AUC) but the conditional interval
excludes zero only for EUR chr1 (+0.0100 [+0.0038,+0.0165]); the chr22 EUR interval is
[−0.0003,+0.0110] (upper bound positive, 4/5 splits) and both 3pop intervals span zero.
`oracle_q` is essentially identical to `oracle_single` (0.842 vs 0.842 for 3pop, 0.993 vs
0.993 for EUR chr22), so fusing the conditioned channels across q does not strengthen
this late-fusion passive bar at this size. The matched TPR@5%FPR gap is positive for EUR
(+0.035 chr22, +0.065 chr1) and negative for 3pop at this matched size (−0.004 chr22,
−0.058 chr1). The matched printed-output result is therefore a positive-mean AUC gap on
every cell, established with a positive interval only on chr1 EUR; the chr22 EUR and both
3pop cells are inconclusive at this sample size, and the matched TPR@5%FPR reading is not
established.

## What these results support
1. **J1**: the accumulated q8 attack beats the metric-specific passive bar over four
   constructions (the held-out best of conditioned AF, the conditioned single release,
   and their label-free and calibrated fusions) in all three cohorts on both
   chromosomes with a positive AUC CI; 5/5 AUC-positive splits (chr22 EAS 4/5). The AUC
   margin is small within a population (+0.007 EUR / +0.004 EAS chr22) and larger across
   populations (+0.066 3pop); the TPR@5%FPR gap against the TPR bar is larger
   (EUR +0.040 chr22 / +0.074 chr1) and positive in 5/5 splits for chr22 EUR/3pop and
   chr1 EUR, 4/5 chr22 EAS, and 3/5 chr1 EAS/3pop; the mean interval excludes zero
   everywhere except chr1 3pop.
2. **J2**: the effect is not the public center/scale (af_pooled 0.69–0.81, far below)
   and not repeated sampling of a fixed input (passive flat; solver-seed repeats gain
   0).
3. **J3**: the interface the attack needs is stated per item and demonstrated on a real
   PLINK 2.0 export; the honest private matrix is never read.
4. **J4**: the randomized PLINK solver changes neither the identity residual nor the
   spanned output, so the attack AUROC is unchanged at this problem size.
5. **J5**: the fixed 1e-9 rank cutoff is precision-blind; a perturbation-derived
   reliability floor on the recovered operand cuts the reported rank from 78 to 30 (6 dp)
   and 11 (4 dp) on chr22 EUR (a heuristic count, not a usable-rank theorem).
6. **J6**: the increment is positive (5/5) at every known-member fraction including 0%,
   with the orientation always the mechanism sign; the labels only control the
   conditioning.
7. **J7**: the honest analysis (lambda-GC, power, recovered-causal count,
   off-causal discovery fraction, effect agreement) is unchanged between the fixed and
   the input-varying pass at this study size.
8. **J9**: the matched printed-output comparison gives a positive-mean AUC gap of the q8
   stack over the passive bar built from the same printed PLINK outputs in every cell
   (+0.0054 EUR chr22, +0.0082 3pop chr22, +0.0100 EUR chr1, +0.0060 3pop chr1); the
   conditional interval excludes zero only for EUR chr1 (+0.0038,+0.0165), and the
   late-fusion passive bar does not strengthen when given the q outputs.

## What these results do NOT support
- **Not a biobank-scale approx claim**: PLINK's randomized `approx` algorithm executes,
  but at n = 240 its output agrees with the exact solver (identical eigenvalues and
  spanned output), so the solver-robustness result is a lower bound on the approx stress.
- **The within-population margin is small** and at the AUC ceiling; the practical
  reading is TPR@5%FPR, where the q8 stack rises from a TPR-bar baseline of 0.960
  (chr22 EUR) / 0.436 (chr22 3pop) / 0.983 (chr22 EAS) / 0.926 (chr1 EUR) / 0.460
  (chr1 3pop) / 0.980 (chr1 EAS) to q8 TPRs of 1.000 / 0.580 / 1.000 (chr22) and
  1.000 / 0.475 / 0.997 (chr1). It is the FLOAT64 cross-population cohort that carries
  the large, interpretable AUC increment; the matched printed-output cross-population
  superiority remains unresolved (J9).
- **q4 is unreliable and sometimes negative in the cross-population and EAS cohorts**;
  the reliable checkpoint is q8.
- **J5 is a reliability heuristic**, not a certified usable-rank theorem; the floor is
  an offline construct from the float64 operand, and on the cross-population arm the
  reliable-subset AUROC drops more (0.806→0.711 at 6 dp), so some removed directions
  carried signal there.
- **J7 is a summary-level stability statement at one study size**; the power values are
  low, the 3pop off-causal discovery fraction (1.9e-4) can include LD proxies of causal
  markers, so it is not a null-error estimate, and there is no per-release cost or
  input-consistency arm.
- **The matched printed-output result is inconclusive in three of four cells**: J9's
  positive-mean AUC gap has a conditional interval excluding zero only for EUR chr1
  (+0.0038,+0.0165); chr22 EUR is [−0.0003,+0.0110] (4/5 splits) and both 3pop cells span
  zero. The matched TPR@5%FPR reading is negative for 3pop at this size.
- **Independent-subject replication (J8) has not been run**; chr1 reuses the same 1000G
  subjects and is not an independent cohort. HGDP availability was not obtained.
- **The interface remains one software package** (PLINK 2.0) and one 1000G release.
