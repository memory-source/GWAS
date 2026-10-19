# EXPERIMENT_RESULTS — Audit-Governed Multi-Site GWAS

Scripts `code/`, raw JSON `results/raw/`, figures `paper/figures/`. Env:
`/data/lvvyongji1/GWAS/env`. Simulated panels + real 1000 Genomes chr22.

> **Round-5 re-run (2026-10-03).** Cross-model review found that the earlier
> scripts used an inconsistent `n_snps` (3000 vs 4000), so the frontier table
> and the mechanism/robustness tables were not on the same panel. All scripts are
> now unified at `n_snps=4000` and re-run; every number in the paper comes from a
> single consistent panel configuration. Real-data validity (S7) added.

## S0 — sanity (`s0_sanity.json`, `s0b_sanity_seeds.json`)
8-seed: uncorrected λ_GC median 4.15 (range 1.17–5.85), PC-corrected median 1.02.
AF-LRT AUROC 0.844; study-fitted PC-distance 0.12 (naive) / ≈0.88 (adversarial);
noise 0.496. Meter valid.

## S1c — frontier, disjoint non-members, 5 seeds, n_snps=4000 (`s1c_frontier.json`)
| tier | AUROC | TPR@1% [95% CI] | structure ctrl |
|---|---|---|---|
| F2_topk | 0.590 | 0.023 [0.009,0.044] | 0.502 |
| F4_round1 | 0.622 | 0.031 [0.014,0.054] | 0.502 |
| F3_round2 | 0.840 | 0.183 [0.119,0.255] | 0.502 |
| F5_full | 0.845 | 0.193 [0.132,0.254] | 0.502 |
| F5_sitewise | 0.964 | 0.588 [0.490,0.696] | 0.502 |

Per-ancestor uniform (F5: 0.860/0.849/0.851).

## S2c — PC-diagnostic leak mechanism, 5 seeds (`s2c_mechanism.json`)
base 0.839; pcdist-only 0.881; pcalign 0.501; base+M2 0.941 (Δ+0.102);
T2 (ref-fit loadings) Δ+0.003; T3 (study-mean center) Δ+0.102; eig ratio 1.000;
permutation 0.486.

## S2b — intermediate reveals, 3 seeds (`s2b_intermediate.json`)
base 0.844; +λ_GC/+cond-bit/+missingness 0.844 (0.000); +PC 0.942 (Δ+0.098).

## S8 — functional-form robustness, 5 seeds (`s8_families.json`)
Same information channel, three functional forms (NOT three channels):
F5: AF-LRT 0.845 / AF-corr 0.833 / ML-dev 0.832 / PC-distance 0.882.
F4(1dp): 0.622 / 0.660 / 0.659 / (PC 0.882) — the LRT is most rounding-fragile.

## S7 — real-data external validity (`s7_real.json`)
1000G Phase-3 chr22, EUR/EAS/AFR pseudo-sites, disjoint 3-way split, MAF
[0.05,0.95], distance-thinned (15 kb), 5 split seeds: AUROC F5_full 0.684,
F3(2dp) 0.681, F4(1dp) 0.564, F2_topk 0.640, F5_perpop 0.808; reference-basis
structure control 0.504. Ordering reproduces at smaller magnitude. Study-fitted
PC-distance on real data is weak and same-signed but attenuated (0.347 naive /
≈0.65 flipped) vs simulated 0.12/0.88; not interpreted as evidence for/against
the PC channel on real data.

## S5 — site-size × precision boundary (`s5_boundary.json`) TPR@1%, n_snps=4000
| n/site | 1 dp | 2 dp | full |
|---|---|---|---|
| 50 | 0.270 | 0.993 | 0.993 |
| 100 | 0.183 | 0.745 | 0.723 |
| 200 | 0.062 | 0.307 | 0.301 |
| 400 | 0.040 | 0.133 | 0.136 |
| 800 | 0.015 | 0.038 | 0.042 |

## S4 — msprime cross-check (`s4_msprime.json`)
full 0.642 / 1dp 0.649 (rounding within noise at that operating point).

## S6 — dimension sweep (`s6_dim.json`)
Δ always positive (+0.006…+0.136); PC-distance beats AF release at all 9 points.

## S3b — sanitized-diagnostics paradox (scripted analyst)
probes 3→0 both families at equal measured quality. Model-level conjecture.

## Defects fixed across review rounds
1. Self-referential negatives → disjoint third non-member panel.
2. maf() dropped AF>0.475 → fixed; tested SNPs ~1874→~3814.
3. Ancestry control was a no-op → reference-basis structure-only control (0.502).
4. Inconsistent n_snps (3000 vs 4000) across scripts → unified at 4000.
5. GT parsing in S7 split only 9 fields → full split + both alleles.

---

# Round 11 — limitation closure (2026-10-06)

Plan: `refine-logs/EXPERIMENT_PLAN.md`; designs: `refine-logs/ROUND11_DESIGN.md`.
CPU env `/data/lvvyongji1/GWAS/env`; GPU `/opt/miniforge/bin/python3`.

## B1 — Per-individual privacy metric (S38) — `results/raw/s38_privacy_metric.json`
5 seeds. pc_loadings AUC .882 / A* .610 / eps_lb(.01) 5.0; af_full .845/.546/4.1;
af_2dp .840/.532/4.2; af_1dp .622/.187/1.1; af_sitewise .964/.805/5.8;
pc_eigonly control .500/.043/0.24; pc_refbasis control .502/.033/0.19.
A*−(2AUC−1)<0 on the leaking channels and >0 on the two at-chance controls
(non-binormal scores). eps_lb is a finite-sample LOWER bound with the false-
negative resolution floored at 1/n_pos (n_pos=n_neg=900; resolvable δ≥~0.01).
For the four high-AUROC channels the `eps_lb_*_valid` flag is False because a
single zero-observed-FN point makes the strict empirical trade-off infeasible
for any finite ε; the 1/n_pos floor is the convention that resolves it, so those
values are read as ordering information.

## B2 — Built sole-egress binding (S39 + enforce.py) — `results/raw/s39_enforcement.json`
library prototype undeclared write: file created (frac 1.00); isolated child
natural-egress: blocked 1.00; sanctioned undeclared push: blocked 1.00; declared
flow 1.00; chain verified 1.00; transported governed af .624 / pc .509.

## B3 — Joint release-bundle leak (S40) — `results/raw/s40_bundle.json`
bundle .962 vs max single .947 (af_site); excess +.015, 5/5 seeds; LOO af_site
.886; af-only .839; af+noise .839 (structural-zero control).

## B4 — Multi-chrom phased IBD + realistic missingness (S41, S42)
S41 (`s41_phased_ibd_genome.json`): chr1/2/22, 50k markers, cM statistic; real
LSR .456/.479/.424 ≤ AF .866/.864/.838; planted 10-cM control .715/.729/.702.
S42 (`s42_missingness_realistic.json`, 5 seeds): rho 0→1 AUROC .490/.629/.755/.893/1.000;
rho=0 null .490; candonly .529/.651/.760/.888/1.000; bonus sweep .49→1.00; super_pop-batch equivalent.

## B5 — Semi-synthetic biobank cohort (S37) — `results/raw/s37_realsynth_scale.json`
Li-Stephens cohorts from real 1000G haplotypes, nsp 150→12000; AF leads PC at
every cell (delta −0.026…−0.088, 5 seeds per cell); in-sample gap →0 with n. Controls: LD-decay,
AF-spectrum, top-10 eig-fraction match the donor library at all cells;
relatedness (parsimony-phase longest match) within ~10% of donor at 6/7 cells
(nsp150 141 vs 128 markers, +11%; at nsp150 cohort CoV 186/141=1.32 vs donor
48/128=0.38, other cells' cohort CoV ~0.3) — the S37 headline does not consume
the IBD channel.

## Integrity audit (EXPERIMENT_AUDIT.md)
2026-10-06, pi-review/glm-flash ultra, read-only. Verdict WARN. Five items,
all addressed: (1) `common.auroc` tie handling fixed to midranks and s41 re-run
(the constant-indicator `auc_cm*` columns no longer store 0.0 artefacts);
(2) s40 sitewise non-member reference index aligned with s38 (conservative
before, now consistent); (3) S43 control reading reworded to the per-control
values; (4) s40 dead imports and the unimplemented ancestry-only control removed
from the docstring; (5) `ROUND11_DESIGN.md` grid text corrected to the delivered
grid. Result-file numbers unchanged: all 30+ quoted round-11 values matched their
JSON keys in the audit.

## Loop-1 (auto-review) fixes
pi-review/glm-flash, round 1, 7.5/10 "almost". Addressed: s38 `eps_lb` recomputed
with an FN resolution floor (1/n_pos) so the value is no longer clip-driven
(pc_loadings 5.0, af_full 4.2, af_1dp 1.1 at δ=0.01) and n_pos=n_neg=900 stated;
paper's A*-vs-2AUC−1 direction corrected (below for leaking channels, above for
controls); §replication's stale "future work" IBD sentence redirected to §scale;
Limitations cohort cross-ref fixed; base-missingness sentence reworded to
"stays saturated".


---

# Round 12 — residual-limitation closure (2026-10-07)

Plan: `refine-logs/EXPERIMENT_PLAN_R12.md`. Closes the four gaps left by the
round-11 result-to-claim judge (`CLAIMS_FROM_RESULTS.md`) and the integrity audit.

## B1 — Per-individual metric on real + semi-synthetic data (S44) — `results/raw/s44_individual_real.json`
The s38 estimators re-run on real 1000G chr22 (8 splits) and an s37 cohort
(nsp=2000, 3 seeds). Real chr22: af_full AUC .727 / A* .350; af_1dp .570/.129;
af_sitewise .865/.584; pc_loadings .665/.271; pc_eigonly control .502/.067;
pc_refbasis control .500/.053. Cohort (3 seeds): af_full .590/.133; af_sitewise
.639/.204; pc_loadings .520/.040; controls ~.496/.017. Ordering preserved on real
(PC<AF, as s16) and on the cohort (AF leads PC); controls at A*≈0 on both. The
f-DP `eps_lb` carries a confidence-corrected `rule3` floor (3/n_pos): valid=True
on the cohort (all channels); on real chr22 the high-AUROC channels still hit a
zero-observed-FN corner, so they are reported as ordering-only there.
Also: `code/s38_privacy_metric.py` now emits both the resolution and rule-3
eps_lb; the simulated values are unchanged for the leaking channels.

## B2 — Biobank-corner cohort cell (S45) — `results/raw/s37_realsynth_scale.json`
One cell added at nsp=16000 (M=3000, 3 seeds; 167 s): af_only .529 vs
pc_calibrated .503, delta −0.026, refbasis .500, eigonly .500, in-sample gap
.003, fst .096. AF still leads and all four controls pass at the largest
feasible cell on a 24 GB GPU (nsp=33000 needs 27 GB for the Gram eigenvectors).

## B3 — Second bundle manifest + uncoarsened single (S46) — `results/raw/s46_bundle_manifest2.json`
Manifest 2 (per-site AF at FULL precision + top-10 eig + lambda_GC@2dp +
call-rate table; pooled AF and p-values omitted): bundle .973 vs coarsened single
.947 (excess +.026) and uncoarsened single .964 (excess +.008, positive 5/5
seeds). Structural-zero control: af_site+noise .964 == af_site .964. LOO
per-site AF drop .233. The excess survives an uncoarsened comparator and a second
manifest, so round-11's +.015 is not a coarsening artefact or a one-manifest
single point (the excess is smaller, +.008, against the uncoarsened single).

## B4 — Native-code egress boundary test (S47 → merged into S39/enforce.py) — `results/raw/s39_enforcement.json`
`code/enforce.py` child now also attempts `ctypes.CDLL(None)` when
`INJECT_NATIVE=1`; the PEP-578 hook denies it on `ctypes.dlopen`/`ctypes.dlsym`.
Result: `native_egress_blocked_frac` 1.0 (3/3 seeds), alongside the unchanged
bypass matrix (library leak 1.0; natural-egress blocked 1.0; sanctioned blocked
1.0; declared flow 1.0; chain verified 1.0; governed af .624 / pc .509). The
stated "no native code" threat-model boundary is now tested, not asserted.

---

# Round 13 — real-data structure + size (2026-10-07)

Plan: `refine-logs/EXPERIMENT_PLAN_REAL.md`. Extends the real-data arm from the
mixed-super-population (AF-first) point to single-population cohorts where the
PC channel dominates, and exposes a study-size confound with a scale-matched rerun.

## B1 — Real-data structure ladder (S48) — `results/raw/s48_real_structure.json`
Real 1000G chr22, identical attack code across cohorts of increasing measured
structure and size (5 split seeds each):
- within a single 1000G population (one group, n_study 28–36): PC dominates AF,
  gap +0.050…+0.134 (GBR .896/.998, CHB .863/.996, JPT .884/.987, TSI .861/.995,
  CDX .884/1.000, YRI .950/1.000, …); 12 single-pop cohorts all gap>0.
- 2-pop mixes (Fst .004–0.020, n_study ~52–66): gap +0.112…+0.172.
- 3/5 super-populations (Fst .075–.087, n_study 450/715): gap −0.005 / −0.030.
Controls: reference-fitted basis 0.51–0.60, eigenvalue-only 0.51–0.60, never
tracking the gap.

## B2 — Real within-population mechanism + granularity (S49) — `results/raw/s49_real_mechanism.json`
Within-population real cohorts (GBR/CHB/YRI on chr1/chr2). On the
PC-dominant real operating point: top-3 loadings recover most of the leak
(0.88–0.96), full loadings 0.99–1.00; eigenvalue-only 0.53–0.60, reference-fitted
basis 0.55–0.60, and an out-of-sample (held-out-half) basis transfers at chance
(0.55–0.58) while the same-half basis is 0.998–1.000 — the in-sample mechanism,
on real data, at the dominant magnitude. (k=10 on n_study 30–36, so k/n≈1/3.)

## B3 — Multi-chromosome within-pop replication (S50) — `results/raw/s50_real_structure_multichrom.json`
chr1/chr2 replicate the within/across split: within-pop gaps +0.051…+0.086 (all
positive); across-population (3/5-pop) gaps −0.039…+0.023 (chr2 3-population is a
mild PC>AF exception at +0.023, consistent with its proximity to the crossing);
controls at
chance on every chromosome.

## B4 — Scale-matched structure ladder (S48, wrap-up) — `results/raw/s48_scalematched.json`
The B1 ladder confounds structure with study size (single-pop n=28-36, 2-pop ~52-66, 3/5-pop 450/715). A scale-matched rerun (`cohort_cell(..., per=20)`, 5 seeds) holds per-group n=20:
- single GBR/CHB (Fst 0): gap +0.091 / +0.076 (PC>AF)
- 3pop (Fst 0.0997, n=60): gap **+0.161** (PC>AF)
- 5pop (Fst 0.0885, n=100): gap **+0.111** (PC>AF)
So at matched small n the PC channel dominates even across populations; the AF-first switch appears only once study size grows (3/5-pop at n=450/715: gap −0.005/−0.030). The real-data ordering is therefore governed jointly by study size and structure — the same small-(n×M) boundary the simulator locates (Section sec:scale) — not by structure alone. The simulator's size-resolved crossing brackets (s21_sensitivity.json, n_snps=8000): [0.027,0.052] @ n=1600, [0.052,0.077] @ n=800, [0.103,0.151] @ n=300 — the crossing value of structure rises as n shrinks.

## Claim check (round 13) — `refine-logs/pi_review_verdict_round13_claimcheck.md`
pi-review/glm-flash ultra: claim_supported YES (high). C7 (within-pop PC>AF across 18 cohort-cells) decisive; C8 (s49 in-sample mechanism, no leak) decisive; C6 supported as a joint size-and-structure gating statement. Text-level refinements applied: quote s21's per-size crossing brackets instead of the single ≈0.04–0.06 interval; state the chr2 3-pop cell as "at the crossing" (chr1 −0.006 there); surface the s50 chr1:YRI eigonly 0.618 (range 0.55–0.62).

## Round 18b — cross-family (Codex/gpt-6-astra) defect closure — 2026-10-08

Codex round C1 (4/10) then C2 (5/10) found engineering defects. Fixed and re-run:

- **S41 estimator (was crashing).** `_cm_run_vs_one` read `cm[len(cm)]` when the final
  marker mismatched (NumPy evaluates `cm[sidx]` before `where`). Clipped `sidx` to
  `m-1`; unit-tested the final-mismatch / all-match / mid-run cases. Full 3-chrom
  re-run (5055 s): real LSR 0.512/0.520/0.471 (chr1/2/22), planted 10-cM control
  0.769/0.766/0.749, AF 0.866/0.864/0.838. JSON regenerated; paper updated.
- **Orientation leak.** `_orient_auroc` selected the sign on `y[h]` but evaluated all
  rows including `h`. Now evaluates only `~h`, and all channels on the same rows share
  one `cal_mask` (paired AF/PC−AF). S48/S49/S50/S48-scalematched re-run. s49 now also
  saves per-field std. New s49: top3 0.907–0.930, full 0.996–0.999, eig 0.415–0.497,
  refbasis 0.416–0.520, oos 0.442–0.601. s48 3/5-pop gaps +0.007/−0.041; scalematched
  3pop +0.193 / 5pop +0.137.
- **S38 language.** Removed "finite-sample lower bound" / "confidence-corrected
  quotable" / ROC-shape inference from the binormal KS.
- **S37.** `combined` is a joint AF+PC logistic; `n_fit`/`fit_frac` recorded; the
  nsp=12000 cell re-run to 3 seeds; new `s37_fitsens.py` (n_fit 5000 vs 20000) shows
  the AF-first gap is fit-size-stable (−0.034/−0.034).
- **S42** described as a generative sensitivity curve (skeleton supplies marker set +
  labels, not measured missingness).
- **enforce.py** child stderr → DEVNULL; s39 re-run all 1.0.

## Round 15 — limitation closure (2026-10-09) — plan: EXPERIMENT_PLAN_LIMITS.md

Addressed the five limitations stated in the 20-page draft, each to a stated
scope. All results below are
raw JSON in `results/raw/`.

- **L1 — real low-structure / fixed-composition size ladder (`s51_real_lowstruct.py`,
  `s51_real_lowstruct_chr{1,2,22}.json`).** Two arms. Pooled within-super-population
  (single group, measured structure 0) at the largest split n: chr22 EUR/EAS/AFR
  gaps +0.075/+0.051/+0.044 (n=167/168/220), chr1 +0.110/+0.165/+0.048, chr2
  +0.181/+0.111/+0.154; both basis
  controls at chance (0.47–0.54). Fixed-composition 3-super-population ladder
  (composition held 1:1:1, only total size changes): chr22 per-group n=40/80/120/160
  gives gap +0.120/+0.042/+0.007/−0.026 at total study size 120/240/360/480 at
  ~constant Fst 0.087–0.091; chr1 +0.058/+0.002/−0.010/−0.018 at Fst≈0.10; chr2
  +0.109/+0.059/+0.018/+0.006 at Fst≈0.086. So at
  fixed superpopulation proportions, jointly increasing the study, reference and
  evaluation panels produces a decreasing split-mean gap on real genomes; the
  crossing locations on chr22/chr1 (n≈260–400) are interpolations of the
  five-split means (top rungs within split noise of zero; chr2 stays positive), so it
  is a joint size effect (the reference panel grows with the study panel), not an
  isolated study-size effect; a zero-structure real cohort at that size is not
  available inside 1000G (single-super-population pools cap at ~500 samples).
- **L2 — eps rule-of-three floor + inversion (`s38_privacy_metric.py`).**
  `validate_tradeoff` now uses the SAME FNR floor as `eps_lower_bound`. The
  rule-of-three floor (3/n_pos) or the default 1/n_pos keeps `eps_lb_rule3_*_valid`
  = True at every delta on every channel; the flag checks consistency with the
  stated floor, not statistical validity. Inversion example
  recorded: af_full (AUC 0.845, eps_lb 4.12) vs af_2dp (AUC 0.840, eps_lb 4.19) —
  higher AUROC, lower eps.
- **L3 — second haplotype/LD estimator (`s52_haplo_second.py`).** Two further
  leave-self-out statistics (max whole-haplotype IBS share; count of released
  haplotypes with IBS share > 0.78). Real unrelated arm at chance (chr22 0.488/0.496,
  chr1 0.497/0.497, chr2 0.497/0.498); a planted 1500-marker copied segment is
  detected (chr22 0.798/0.797, chr1 0.776/0.788, chr2 0.772/0.779). The haplotype
  negative is not specific to the
  longest-shared-run estimator.
- **L4 — fixed eigenvalue-informed basis (`s53_eigbasis.py`).** Fixed,
  data-independent bases (identity, cosine, and a λ-weighted basis that reuses the
  released spectrum) are all at chance on the simulated and real within-population
  cells (0.46–0.53) on chr22/chr1/chr2 × EUR/EAS/AFR, while the study-fitted loadings recover the leak (0.79–0.90). The
  eigenvalue-only negative does not rest on the basis being random.
- **L5 — fit-fraction sweep (`s37_fitsens.py`, 5 seeds).** At nsp=8000 the AF-first gap is
  −0.034 (sd 0.008) at every n_fit in {5000,10000,20000,24000 (full)}; at nsp=16000
  it is −0.026 at n_fit ∈ {5000,10000,20000} (the full 48000-row fit
  was not run). The
  subsampled-fit cells are representative.

---

# Round 15b — real-gap uncertainty + density + seeds (2026-10-11/13)

## S54 — real-gap permutation null + per-split CI (`s54_real_permnull.py`)
Recomputed AF/PC on the s48 cohorts with the identical attack code and added a
stratified label-permutation null (5000 permutations on the held-out evaluation
half) plus a bootstrap CI over the five per-split gaps. `s54_real_permnull.json`:
- within:GBR gap +0.085 (sd 0.041), CI [+0.058,+0.126], p_gt0 1.00, pc_p 2e-4;
- within:CHB +0.108 (sd 0.093), CI [+0.047,+0.202], p_gt0 1.00, af_p 1.8e-2;
- within:YRI +0.053 (sd 0.018), CI [+0.038,+0.070], p_gt0 1.00;
- 3pop +0.007 (sd 0.016), CI [−0.008,+0.020], p_gt0 0.60 (near parity, straddles zero);
- 5pop −0.041 (sd 0.042), CI [−0.077,−0.005], p_gt0 0.20 (AF-first).
Observed AUROCs reproduce the s48 JSONs to <0.002 (same attack code), so the
dispersion is the split-to-split variability of the same statistic.

## S55 — real within-population density axis (`s55_real_density.py`)
GBR/CHB/YRI chr22 at three thinning tiers (`s55_real_density.json`):
- D_thin (~17.3k markers): gaps +0.007 / +0.035 / +0.000
- D_med (~9.8k markers): +0.005 / +0.038 / +0.000
- D_dense (~46.9k markers): +0.002 / +0.021 / +0.000
PC >= AF at every tier; the within-population ordering is density-stable (both
basis controls at chance 0.43–0.58). AF alone saturates toward 1.0 in these tiny
cohorts, so the gaps are small absolute differences.

## S37b / S42b — seed bump 3→5
- s37 (5 seeds/cell): nsp150→16000 gap −0.088/−0.070/−0.064/−0.034/−0.038/−0.031/−0.026
  (sd up to 0.027); fitsens nsp8000 −0.034 (sd 0.008) flat across n_fit, nsp16000 −0.026.
- s42 (5 seeds): rho 0→1 AUROC .490/.629/.755/.893/1.000; candonly .529/.651/.760/.888/1.000.

## S56 — real structure-zero switch, fixed-reference design (`s56_real_poolmax.py`)
The pooled single-super-population panel (one group, measured structure zero by
construction) with the STUDY grown and the REFERENCE held fixed at
$n_{\text{ref}}{=}150$ (non-members absorb the remainder; where the pool cannot
supply 150 non-members the reference is capped and its size is recorded). chr22+
chr1+chr2, 5 split seeds (`s56_real_poolmax_chr{22,1,2}.json`); the PC-AF gap
stays POSITIVE at every study size (EUR/EAS/AFR):
- n167 (all pops, n_ref 150): chr22 +0.080/+0.102/+0.095; chr1 +0.140/+0.141/+0.106; chr2 +0.170/+0.116/+0.179
- n350 (AFR n_ref=150): chr22 +0.102; chr1 +0.140; chr2 +0.181 (sd 0.018–0.037)
  (EUR/EAS at n350 have n_ref 147/148 and only 6 non-members; their gaps are positive but high-variance, sd up to 0.41)
Reference-fitted-basis and eigenvalue-only controls stay near chance (0.37–0.63)
in every cell. So the structure-zero PC>AF ordering does not flip inside the
largest real cohorts the single-super-population pool (EUR 503 / EAS 504 / AFR
661) permits with a full-size reference; a zero-structure real AF-first cohort
remains unobserved inside 1000G.

---

# Round 17 — cross-model defect closure after the round-16 audit (2026-10-14)

## S59 — GPU AUROC tie-handling (`code/s59_auroc_tiecheck.py`) — `results/raw/s59_auroc_tiecheck.json`
`common_gpu.auroc_gpu` builds ranks via `argsort` and `(ranks*y).sum()`, so it
assigns arbitrary distinct ranks to tied scores. Constant scores `[1,1,1,1]` with
labels `[1,1,0,0]`: CPU midrank 0.5, GPU rule 0.0. Half-tied 0/1 vector: 0.505 vs
0.252 (gap 0.253). A regenerated archived s28-style cell (n=3000 AF-LRT scores)
has tie_mass 0.0 (3000 unique of 3000), so the archived headline precision is
unaffected; the kernel is disclosed as not tie-correct (not rewritten, to avoid
invalidating the s27-s36 artifacts).

## S60 — monotone share + per-seed crossing (`code/s60_share_crossing.py`) — `results/raw/s60_share_crossing.json`
On the s51 fixed-composition 3-super-population ladder (per-sample SUPER-POP
label), the cumulative squared-projection share is non-decreasing in k on every
rung (per40/80/120/160: e.g. 0.109→0.291 at per40), removing any concern that the
-0.026 top-rung gap is an artefact of a per-component share that rises by
construction. Per-seed interpolation of the chr22 crossing gives seed values
{0:[239,242,412],1:[360],2:[280],3:[383],4:[]} — 3 single-crossing seeds (mean
n≈341, sd≈44), one seed with three crossings and one non-crossing. The reported
crossing bracket (n_study≈260–400) is stated as a five-split-mean interpolation.

## S56 — fixed-reference rewrite (`code/s56_real_poolmax.py`)
The earlier "study-heavy" split let the reference shrink with the study, which
mechanically lowers AF. Rewritten to hold n_ref at 150 (cap only where the pool
cannot supply 150 non-members) and re-run on chr22/chr1/chr2. AFR keeps n_ref=150
through n_study=350 with gap +0.102/+0.140/+0.181 (sd 0.018–0.037); EUR/EAS at
n350 have n_ref 147/148 (with 6 non-members) and positive but high-variance gaps (sd up to 0.41).
Controls near chance (0.37–0.63). The gap stays positive at every study size.

## Text corrections
- s54 scope: the permutation null covers the three within-population cohorts plus
  the 3-/5-population cells (five tested cohorts), and tests each channel's
  separation, not the PC−AF difference; `gap_p_gt0` is the positive-split
  fraction.
- s37 controls: LD decay matches at the five cells whose panel supports the
  shortest bin (undefined, not mismatched, at M=3000/6000); AF spectrum and
  eigenvalue fraction match at all cells; relatedness within ~10% at 6/7 cells.
- s56 wording: "necessarily falls to 30" dropped; the largest full-reference cell
  is AFR at n350; control range 0.37–0.63.
- s42: the two missingness curves differ by ≤0.01 at interior rho (not
  "identical").
- CLAIMS/TRACKER: s56 entry, s37 cell/seed counts, and the stale "permutation
  p-values pending" note corrected.
