# Experiment Plan (GPU / scale round)

**Problem**: The diagnostics-leakage paper's magnitudes all come from small
panels (n<=1600/group, M<=20k), and its real-data arm is capped at 1000G
chr1/chr2/chr22 (~2500 samples, ~5k markers/chr after thinning). Reviewers want
the ordering's regime boundary quantified, the statistics tightened, and the
real-data caveats (marker density, super-population coverage, the deferred
phased-IBD arm) closed with bigger data and more compute.
**Method Thesis**: The released study-fitted PCA spectrum is the dominant
*canonical* intermediate leakage channel; its mechanism is in-sample-ness of the
fitted object; the ordering is governed by measured population structure, and a
correctly-scoped reveal manifest can govern it.
**Date**: 2026-10-05

Runtime facts (measured on the box): 8x RTX 4090 D 24GB, CUDA 12.1, system
python3 torch 2.2.0+cu121 (project conda env is torch-free). Decompression
895 MB/s vs parse ~4.7k variants/s → the VCF bottleneck is the *parse* (CPU),
not a GPU job. Exact `torch.linalg.eigh` Gram PCA fits to n~3e4 (~20GB); n>=6e4
needs randomized SVD (measurably ~0.02 biased LOW → validation-gated).

## Claim Map
| Claim | Why It Matters | Minimum Convincing Evidence | Linked Blocks |
|-------|-----------------|-----------------------------|---------------|
| C1 The PC-vs-AF ordering has a located regime boundary | Converts "small-site phenomenon" from a hedge into a measured scaling law | PC-AF gap over n x M grid at biobank n, fitted zero-contour, plus a scale test of the Fst crossover | B1, B2 |
| C2 Every headline number is tightly significant | Removes the "5-30 seeds / wide CI" objection | 1e5-permutation official nulls + bootstrap CIs; 500-1000-seed paired negation sweep | B3 |
| C3 Loadings (not eigenvalues) carry the leak, monotonically | The governance recommendation depends on the exact k-ladder | Exact k=1..10 + eigenvalue-only arm, many seeds, sim/real/msprime | B3 |
| C4 Real-data conclusions are density- and ancestry-robust | Retires "~5k markers/chr" and "3 super-pops" caveats | Dense-marker tiers on 3 chromosomes; all 5 1000G super-populations | B4 |
| C5 A real phased IBD attack is weak among unrelated individuals | Closes the explicitly-stated phased-arm gap; bounds the haplotype channel | Real phased H0/H1 long-segment attack + planted positive control | B5 |
| Anti-claim: a GPU VCF parser/decompressor is needed for the real arm | Avoids dressing the limitation with a GPU that changes nothing real | Measured parser benchmark (gzip vs pigz vs numba) | B6 |

## Paper Storyline
- Main paper must prove: C1 (located boundary + scale-invariant crossover),
  C3 (exact ladder), C4 (density + 5-pop).
- Appendix can support: C2 (permutation/bootstrap error bars), C5 (phased IBD
  bound), B6 (parser benchmark as a negative result), B7 (nonlinear meta-attack).
- Experiments intentionally cut: semi-synthetic 50k HAPGEN cohort (would need a
  simulator label; does not change a real-data claim); gnomAD (aggregate
  frequencies only — unusable for per-individual scoring); CUDA jax / cupy /
  nvcomp GPU decompression (measured non-bottleneck).

## Experiment Blocks

### Block 1 (B1): GPU biobank-scale boundary grid
- Claim tested: C1.
- Why it exists: the paper cannot see where the PC-dominance boundary is at
  biobank dimensions, nor whether the Fst crossover moves with scale.
- Dataset / split / task: parametric 3-pop panel; disjoint 3-way split;
  n/group in {1e4, 3e4}, M in {1e5, 3e5, 1e6}, Fst in {.02,.04,.08,.12},
  8 seeds/cell.
- Compared systems: AF-only vs train-half-calibrated PC vs logistic-combined.
- Metrics: AF/PC/combined AUROC; paired PC-AF gap + std; in-sample mechanism
  (own vs other basis); fst_like; top10 var-frac; reference-basis control;
  within-ancestry permutation control.
- Setup details: exact blocked-Gram PCA (`torch.linalg.eigh`), TF32 OFF, system
  python3 torch, GPUs 2/3.
- Success criterion: a monotone/contoured boundary; overlap cells reproduce
  s12/s17/s21 within +/-0.005; controls in [0.48,0.52] / [0.47,0.51].
- Failure interpretation: if the boundary is AF-saturation-driven (iso-curve of
  AF ~0.86-0.9), the paper states the frontier precisely and reports PC
  dominance as small-(n*M) only.
- Table / figure target: Fig 8(a) heatmap + zero-contour; a boundary sentence in
  the Discussion.
- Priority: MUST-RUN.

### Block 2 (B2): Fst crossover at scale
- Claim tested: C1 (second half).
- Why: the crossover (measured ~0.04-0.06) was measured only at n<=1600.
- Dataset: reuses B1's Fst axis (near-zero marginal cost).
- Metrics: crossover measured-structure location vs (n, M), bootstrap CI over
  seeds.
- Success criterion: report whether the crossover structure level is
  scale-invariant; either answer refines the scoping.
- Table / figure target: Fig 8(b) crossover-vs-n by M.
- Priority: MUST-RUN.

### Block 3 (B3): GPU permutation/bootstrap + negation + ladder
- Claim tested: C2 and C3.
- Why: the frontier CIs are 300-draw/5-seed; the crossover cells are within
  noise of a sign flip; the k-ladder jumps 1,2,3,5.
- Dataset: frontier tiers (sim); the s12 far-corner cell; the s21 crossover
  grid; real chr22/chr1/chr2 + msprime for the ladder.
- Metrics: 1e5 stratified-permutation nulls (per tier, per ancestor) + bootstrap
  CIs; 500-1000-seed paired PC-AF delta with P(delta>0); k=1..10 + eig-only with
  CIs.
- Setup: batched matmul-form attacks + rank-sum batched AUROC (one permutation
  sweep = one tensor op).
- Success criterion: CIs tighten ~2.5x; the crossover sign is pinned or shown
  flat-zero across the saturated-AF corner; the ladder saturates at a measured k.
- Failure interpretation: a far-corner sign flip ⇒ the paper says "ties AF only
  where AF has saturated".
- Table / figure target: Table 1 CIs; Fig 8(c) ladder.
- Priority: MUST-RUN.

### Block 4 (B4): Real density + 5 super-populations
- Claim tested: C4.
- Why: retires "~5k markers/chr" and "3 super-population pseudo-sites".
- Dataset: chr22/chr1/chr2 at density tiers D5k / D50k; all five 1000G
  super-populations (adds SAS 489, AMR 347).
- Compared systems: s16-style orientation + granularity + per-pop release at
  each tier / pop set.
- Metrics: PC-AF gap and IBS AUROC vs marker density; per-superpop AF/PC AUROC;
  5-pop vs 3-pop paired comparison.
- Setup: compiled (numba/pigz) parser; cache per-chr int8 npy. CPU only (PCA at
  n<=2504 is ms).
- Success criterion: ordering stable across densities; 5-pop does not reverse it.
- Failure interpretation: a density-driven reversal is a real finding; report it.
- Table / figure target: real-data paragraph + appendix table.
- Priority: MUST-RUN.

### Block 5 (B5): Real phased IBD attack
- Claim tested: C5.
- Why: the paper explicitly defers "arm has not yet exploited H0/H1 separately".
- Dataset: dense phased 1000G chr1-3 (H0, H1), 2-bit packed.
- Compared systems: dosage IBS window-sharing (s26) vs phased longest shared run
  (H0/H1) at length thresholds L in {20,50,100}; planted-IBD positive control as
  invariant.
- Metrics: phased IBD AUROC vs AF baseline; null-unplanted; planted control.
- Success criterion: planted control high; real phased IBD <= AF (confirms the
  haplotype channel is weak among unrelated samples) OR finds structure the
  dosage proxy missed.
- Table / figure target: Fig 8(d) / appendix.
- Priority: MUST-RUN (closes a stated gap); GPU kernel at m>=1e5.

### Block 6 (B6): VCF parser benchmark (anti-dressing control)
- Claim tested: the anti-claim.
- Why: turn "we could use a GPU VCF parser" into a measured negative result.
- Compared systems: single-thread gzip loop vs pigz -dc -p16 + parser vs
  numba-vectorized parser vs a GPU string-parsing straw man.
- Metrics: wall-clock and a correctness checksum per path.
- Success criterion: pigz/numba give 20-50x; the GPU path is not faster.
- Priority: NICE-TO-HAVE.

### Block 7 (B7): GPU nonlinear meta-attack
- Claim tested: composition robustness (supporting).
- Why: confirms the 0.000 rows are structural (lambda_GC / condbit are per-panel
  constants) and the IBS feature stays at chance unplanted at scale.
- Metrics: MLP-CV AUROC per feature set at 500 seeds vs the logistic baseline.
- Success criterion: reproduces +0.098 and the three structural zeros.
- Priority: NICE-TO-HAVE.

## Run Order and Milestones
| Milestone | Goal | Runs | Decision Gate | Cost | Risk |
|-----------|------|------|---------------|------|------|
| M0 | GPU harness sanity | common_gpu smoke vs a CPU s12 cell | reproduce within +/-0.005 | 5 min | torch/CPU mismatch |
| M1 | Boundary grid | B1 | controls clean, no OOM | 1-2 GPU-h | n=3e4 eigh OOM |
| M2 | Statistics | B3 | CIs tighten, sign pinned | 30 min | permutation kernel bug |
| M3 | Real data | B4 | density-stable ordering | <30 min | VCF parse |
| M4 | Phased IBD | B5 | planted control high | 1-2 h | dense phased memory |
| M5 | Polish | B6, B7 | negative result recorded | 20 min | none |

## Compute and Data Budget
- Total estimated GPU-hours: ~4-6 GPU-hours, sharded over GPUs 2/3 (+0/1).
- Data preparation: chr1/2/22 local; optional chr1-3 dense phased parse
  (one-time, ~30 min).
- Human evaluation: none.
- Biggest bottleneck: `torch.linalg.eigh` VRAM at n=3e4 (~20GB peak); n>=6e4
  needs validation-gated randomized SVD.

## Risks and Mitigations
- Randomized-SVD bias (~0.02 low on the PC channel): keep exact Gram default;
  validate exact-vs-rSVD at n=3e4; exclude failing cells.
- TF32 silently corrupting fitted loadings: explicitly disabled in common_gpu.
- Real arm cannot exceed ~2500 individuals: state it; do not oversell GPU; the
  semi-synthetic cohort is optional and must be labeled.
- GPU decompression is a non-bottleneck: documented as a negative result (B6).

## Final Checklist
- [x] Main paper tables are covered (boundary, CIs, ladder, real density/5-pop)
- [x] Novelty is isolated (governance prototype + in-sample mechanism)
- [x] Simplicity is defended (eigenvalue-only control, structural-zero rows)
- [x] Frontier contribution justified or explicitly not claimed (GPU used only
      where it helps; parser/decompress GPU explicitly not claimed)
- [x] Nice-to-have runs separated from must-run runs

## Outcome (2026-10-06) — see refine-logs/GPU_ROUND_RESULTS.md
- **B1/B2 (S27, S35):** at biobank scale the AF channel saturates first and PC
  loses everywhere past nsp~500; PC only narrowly leads at very small marker
  counts (M=5k, nsp=500, +0.008). The ref-basis control collapses PC to ~0.50
  and the in-sample gap shrinks monotonically with study size
  (0.51 -> 0.08) — the mechanism is in-sample-ness, not the eigenvalues.
- **B3 (S28, S29):** both channels are far above the stratified permutation null
  (p<1e-4); the small-n PC-AF gap is +0.09 with P(delta>0)=1.0 but decays to
  ~+0.009 at M=8000; the k-ladder saturates near k=8-10 and an eigenvalue-only
  basis stays at 0.51 on sim/real/msprime.
- **B4 (S30, S31):** the real ordering is density-stable (PC stays ~0.66 while
  AF rises with density; gap widens) and 5-pop does not reverse it (chr22
  5pop .708/.611; chr1+2 5pop .780/.611).
- **B5 (S36):** corrected phased IBD: planted control validated (0.671), real
  unrelated-member attack at chance (0.486) — the haplotype/IBD channel is weak
  among unrelated 1000G samples, bounding it below the AF release.
- **B6 (S34):** VCF parsing is parse-bound (pigz pipe 2472 vs gzip 3290
  variants/s); GPU decompression is not a win — negative result recorded.
- **B7 (S33):** an MLP reproduces the logistic result exactly; lambda_GC,
  cond-bit, and unplanted-IBS add 0.000 — the structural zeros are confirmed.
