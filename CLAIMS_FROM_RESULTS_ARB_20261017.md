verdict: yes
confidence: medium
integrity_status: pending (round-22 codex integrity re-audit in progress, thread 01a14948)
round: round 22 — fair comparison and interface confirmation (s67/s68/s69/s70/s71/s72)

Center claim (NARROWED): a cooperating analysis partner that keeps every per-SNP
genotype count unchanged, and makes only bounded input modifications, obtains membership
information beyond the strongest of four constructed passive adversaries (the AUROC
oracle over conditioned AF, the conditioned single release, and their label-free and
calibration-weighted fusions) in the final PCA outputs — established in float64, and
carried to real PLINK 2.0 printed output as a positive-mean AUC gap in every cell of the
matched comparison (s72) whose conditional interval excludes zero only on chr1 EUR. The
approximate-solver, precision, known-member-fraction and scientific-cost findings are
established SEPARATELY. The matched printed-output TPR@5%FPR reading is NOT established.
The larger float64 increment is cross-population; matched printed-output evidence is
positive with an interval excluding zero only for chr1 EUR.

Judge: codex / gpt-6-astra (read-only, cross-family vs deepseek-flash executor),
thread 01a14948-947a-7f22-9862-5dd2d9935d3c. Evidence:
results/raw/s67_joint_passive_chr{22,1}.json, s68_interface_chr22.json,
s69_reliable_dirs_chr{22,1}.json, s70_known_fraction_chr{22,1}.json,
s71_protection_cost_chr{22,1}.json, s72_printed_matched_chr{22,1}.json;
code/s67_joint_passive.py, code/s68_interface.py, code/s69_reliable_dirs.py,
code/s70_known_fraction.py, code/s71_protection_cost.py, code/s72_printed_matched.py;
refine-logs/EXPERIMENT_PLAN_ARB_20261017.md,
refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md.

## Claim map (round-22; extends the deep round's D1-D9)

| Claim | Verdict |
|-------|---------|
| J1 — the accumulated q8 increment survives the metric-specific AUROC/TPR oracle over four passive constructions, in float64 | yes (chr22 AUC +0.0066 EUR / +0.0656 3pop / +0.0037 EAS; chr1 +0.0130/+0.0443/+0.0048; positive AUC CI; 5/5 AUC splits except chr22 EAS 4/5; TPR@5%FPR gap against the TPR-best bar +0.040 EUR / +0.144 3pop / +0.017 EAS chr22 and +0.074/+0.015/+0.017 chr1; 5/5 within-population TPR splits except chr22 EAS 4/5 and chr1 EAS 3/5) |
| J2 — the increment is not the public center/scale and not repeated sampling of a fixed input | yes (af_pooled 0.69–0.81; passive flat; solver-seed same-input span move 6e-8 mean, TPR gain 0.000) |
| J3 — the interface items are demonstrated on a real export (local software capability; re-run/centre-scale access is a threat-model assumption) | yes (s68 manifest: loadings/eigenvalues backed by real PLINK artifacts; center-scale/re-run assumed; honest matrix not read) |
| J4 — the attack survives the approximate (randomized) solver | yes at this problem size (identical eigenvalues and span; identity residual 6.242e-6 both; q8 AUROC 0.998637 both); scoped as a lower bound (n = 240) |
| J5 — the reliable-direction count as a function of export precision | yes as a heuristic count (s69 perturbation-derived floor on the recovered operand: chr22 EUR geometric rank 78 → reliable 30 at 6 dp, 11 at 4 dp; 3pop reliable AUROC drops 0.806→0.711 at 6 dp); the floor is an offline construct from the float64 operand, so this is a descriptive count, not a certified usable-rank theorem |
| J6 — the increment holds across the known-member fraction 0/5/10/50% (fraction of honest members, drawn from the calibration half) | yes (positive 5/5 at every fraction; orientation is the mechanism sign, labels only control the conditioning; EUR +0.064 at 0% down to +0.011 at 50%; passive bar is the reachable label-free z-fusion) |
| J7 — the honest GWAS summary statistics are stable under an input-varying pass | partial (lambda-GC, power, recovered-causal count, off-causal discovery fraction and effect agreement unchanged between fixed and varying pass; power low at this n; 3pop off-causal discovery fraction 1.9e-4, causal proxies can contribute so it is not a null-error estimate; no per-release cost arm) |
| J8 — independent-subject replication | not run (chr1 reuses 1000G subjects; HGDP not obtained) |
| J9 — the matched printed-output comparison (accumulated attack vs the passive bar built from the same printed outputs) | partial (s72: positive-mean AUC gap of q8 over the printed passive bar in every cell — +0.0054 EUR chr22, +0.0082 3pop chr22, +0.0100 EUR chr1, +0.0060 3pop chr1; conditional interval excludes zero only for EUR chr1 [+0.0038,+0.0165], chr22 EUR [−0.0003,+0.0110] 4/5, both 3pop span zero; late-fusion passive bar does not strengthen with all q outputs); the matched TPR@5%FPR reading is not established |

## What the results support
- J1 in all three cohorts on both chromosomes against the AUROC oracle over four
  passive constructions, in float64; the margin is small within a population and larger
  across populations.
- J2: the effect is neither the released center/scale nor fixed-input resampling.
- J3/J4: the interface is stated per item and the attack survives the randomized
  solver at this size.
- J5: the precision boundary — a perturbation-derived floor on the recovered operand
  removes spurious directions (a descriptive count under one plausible floor).
- J6: the calibrated-attacker advantage is not load-bearing; the increment is largest
  when no member labels are used for orientation.
- J9: in the matched printed-output cell the q8 stack gives a positive-mean AUC gap over
  the passive bar built from the same printed PLINK outputs in every cell, with a
  conditional interval excluding zero only on chr1 EUR; the late-fusion passive bar does
  not strengthen when it is given the q outputs.

## What the results do not support
- The matched printed-output TPR@5%FPR reading (J9): the cross-population matched TPR gap
  is negative at this size.
- The matched printed-output AUC gap is supported as a positive mean in every cell but is
  established with a positive interval only on chr1 EUR; the chr22 EUR and both 3pop cells
  span zero at this sample size.
- Not a biobank-scale approx claim: PLINK's randomized `approx` executes but at n=240
  agrees with the exact solver, so the result is a lower bound.
- The within-population margin is small and at the AUC ceiling; the reliable checkpoint
  is q8.
- J5 is a reliability heuristic, not a certified usable-rank theorem.
- J7 is a summary-level stability statement at one study size.
- J8 (independent subjects) not run.
- One software package, one 1000G release.
