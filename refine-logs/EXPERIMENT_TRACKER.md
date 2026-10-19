# EXPERIMENT_TRACKER — Audit-Governed Multi-Site GWAS

Env: `/data/lvvyongji1/GWAS/env`. All scripts unified at n_snps=4000 (round-5).

| Run | Script | Status | Key result |
|---|---|---|---|
| S0  | s0_sanity.py | DONE | λ_GC uncorr 1.17→corr 0.95; attack 0.836 |
| S0b | s0b_sanity_seeds.py | DONE | λ_GC uncorr median 4.15 / corr 1.02 (8 seeds) |
| S1c | s1c_frontier.py | DONE | F2 0.590 / F4 0.622 / F3 0.840 / F5 0.845 / per-site 0.964; ctrl 0.502 |
| S2c | s2c_mechanism.py | DONE | base 0.839 → +PC 0.941 (Δ+0.102); T2 Δ+0.003 |
| S2b | s2b_intermediate.py | DONE | base 0.844 → +PC 0.942 (Δ+0.098) |
| S8  | s8_two_family.py | DONE | AF-LRT 0.845 / corr 0.833 / ML 0.832 / PC 0.882 |
| S7  | s7_real.py | DONE | real 1000G chr22, 5 split seeds: F5 0.684 / 2dp 0.681 / 1dp 0.564 / top-200 0.640 / perpop 0.808; ctrl 0.504 |
| S5  | s5_boundary.py | DONE | n50 2dp 0.993 → n800 1dp 0.015 |
| S6  | s6_dim.py | DONE | Δ always positive (+0.006…+0.136) |
| S3b | s3b_paradox.py | DONE | probes 3→0 both families, quality parity |
| S4  | s4_msprime.py | DONE | full 0.642 / 1dp 0.649 |
| S9  | s9_orientation.py | DONE | real PC naive 0.652 / flipped 0.348, calibrated 0.652 < AF 0.688 (below AF, 20/20 splits); sim calibrated 0.910 > AF 0.794; orientation 'farther' on real (farther_frac 1.0) |
| S10 | s10_paired.py | DONE | 25 seeds, held-out orientation: PC 0.879 vs AF 0.832, paired +0.046 [0.037,0.055] p<1e-4; combined +0.103; sign 25/25 |
| S11 | s11_robust.py | DONE | S9d in-sample 0.88 / out-of-sample 0.50 / swapped-null 0.51; S9c eig-only 0.51, full 0.88; S9b ctrl 0.502→0.63 sensitivity; S9a batch-missingness 0.76 |
| S12 | s12_crossover.py | DONE | PC−AF increment positive across small/mid panels, peaks +0.15 (400,1000), crossover −0.006 at (1600,8000) |
| S13 | s13_real_loadings.py | DONE | real granularity: k1 0.51 k3 0.55 k5 0.58 full 0.65, eig-only 0.51 (loadings needed on real too) |
| S14 | s14_real_mechanism.py | DONE | real held-out-site: own 0.685 vs other 0.516 / swapped 0.516 |
| S15 | s15_controls.py | DONE | batch-missingness coverage controls (full 0.846 / uniform-70 0.800 / batch 0.759); Fst-like + top10 var-frac |
| S16 | s16_multichrom.py | DONE | chr1/chr2 replication: PC < AF both (0.675/0.689 vs 0.768/0.752); mechanism + granularity replicate |
| S17 | s17_scale.py | DONE | site-size sweep M=10k: n1200 AF 0.909/PC 0.897, n3000 AF 0.807/PC 0.745; +PC incr stable +0.056 |
| S18 | s18_msprime_pc.py | DONE | msprime PC cross-check: PC 0.691 > AF 0.627; granularity ladder + in-sample mechanism replicate |
| S19 | s19_haplotype.py | DONE | haplotype/LD channel 0.990 >> AF 0.845; OOS transfer 0.487 (in-sample-only) |
| S20 | s20_extra.py | DONE | M3 missingness channel graded 0.49→1.00 vs call-rate excess; null 0.506; structure comparison |

Cross-model code-review defects fixed: self-referential negatives; maf() AF>0.475
drop; no-op ancestry control; inconsistent n_snps; S7 GT parsing.

Round-8 (result-to-claim gap closure) new scripts: real_panel.py (cached parse),
s9/s10/s11/s12/s13. Key outcome: PC dominance is simulation-scoped; real data
reproduces the PC channel at sub-AF magnitude with the same loadings-granularity
mechanism; PC-vs-AF paired gap now formally significant; batch-missingness
channel quantified; crossover map locates PC−AF ≈ 0 only at AF saturation.

Round-9 (limitation-closure / bigger data) new scripts: real_multi.py (chr1/chr2
loader), s16/s17/s18/s19/s20. Key outcome: real ordering replicates on chr1 and
chr2; msprime confirms PC dominance and mechanism; a distinct haplotype/LD
channel (0.99) is also in-sample-only; M3 missingness is a real graded channel
when membership-correlated; PC dominance does NOT hold at realistic site sizes
(n≥1200) where AF saturates — the honest boundary of the headline claim.

Paper: `paper/main.tex` / `main.pdf` (6 pages); compiles with 0 errors (a few
overfull hboxes remain).

Round-10 (limitation-closure / governance + model-dependence) new scripts:
governance.py (manifest + egress monitor + hash-chained ledger), real_phased.py
(keeps VCF phase), s21_sensitivity, s22_admixture, s23_governed, s24_crossrun,
s25_paradox, s26_ibd.

| Run | Script | Status | Key result |
|---|---|---|---|
| S21 | s21_sensitivity.py | DONE | PC−AF gap decreases MONOTONE in Fst (all 3 n): +0.05 (Fst .02) → −0.08 (Fst .24); crosses 0 at Fst≈0.10 = real-1000G structure |
| S22 | s22_admixture.py | DONE | PC dominance needs discrete structure: alpha 0.2 → PC−AF −0.115; 10 → +0.096; msprime admix +0.13 |
| S23 | governance.py + s23_governed.py | DONE | manifest/egress/ledger real; 5 injections blocked + 1 coarsened, chain verified; governed AF 0.622 / PC 0.51 vs ungoverned 0.845 / 0.882 |
| S24 | s24_crossrun.py | DONE | cumulative composition 0.941 > max single 0.882 (5/5); anti-composition clause flags the refinement |
| S25 | s25_paradox.py | DONE | paradox upgraded to result: dE +0.24 at p_san=1, 0 when p_san>=p_stop; null_fixed control shows NO paradox (falsifiable) |
| S26 | real_phased.py + s26_ibd.py | DONE | phased/IBS attack: planted 0.866, unplanted 0.486, real chr22 0.463 (<= AF 0.653); 0.99 adjacent-pair = fitted-table artifact |

Note: 1000G VCFs are PHASED (GT uses '|'); real_panel.py/real_multi.py discard
the phase — real_phased.py keeps it for the IBD attack.

Paper: `paper/main.tex` / `main.pdf` (12 pages); compiles with 0 errors.

## GPU / scale round (2026-10-05) — plan: EXPERIMENT_PLAN_GPU.md

Env: system python3 (torch 2.2.0+cu121) at /opt/miniforge/bin/python3 for GPU;
project env for CPU. GPUs 0-5, sharded. Panels HOST-resident int8, PCA/attacks
stream marker blocks to GPU (peak GPU = O(n^2) Gram).

| Run | Script | Purpose | Priority | Status | Notes |
|---|---|---|---|---|---|
| G0  | code/common_gpu.py | GPU kernel harness (sim/AF-LRT/PC/exact-Gram PCA/rank AUROC) | MUST | DONE | smoke vs CPU within 1e-5 |
| S27 | code/s27_gpu_scale.py | biobank boundary grid (nsp x M x Fst), exact k-ladder, ref-basis + in-sample controls | MUST | DONE | boundary: AF dominates past nsp~500; in-sample gap collapses 0.51->0.08 with n |
| S28 | code/s28_permnull.py | stratified-permutation nulls (all p~5e-5) + bootstrap CIs | MUST | DONE | both channels highly significant; gap CI tight |
| S29 | code/s29_negation.py | 600-seed paired delta + exact k=1..10 ladder (sim/real/msprime) | MUST | DONE | sim delta +0.09 P=1.0; ladder saturates k~8-10; eigen-only 0.51 |
| S30 | code/s30_density.py | real density tiers (chr22/1/2) | MUST | DONE | ordering stable; PC-AF gap grows with density (-0.06->-0.26 chr1) |
| S31 | code/s31_fivepop.py | all 5 1000G super-pops (chr22+chr1/2, 3pop vs 5pop matched) | MUST | DONE | 5-pop does not reverse: chr22 5pop .708/.611, chr12 5pop .780/.611 |
| S32 | code/s32_phased_ibd.py | phased H0/H1 long-segment IBD | MUST | SUPERSEDED | self-match design flaw -> s36 |
| S33 | code/s33_meta_attack.py | GPU nonlinear meta-attack (structural-zero check) | NICE | DONE | combined 0.985 vs PC 0.799; MLP==logistic; lambda/condbit/IBS add 0.000 |
| S34 | code/s34_parse_bench.py | VCF parser benchmark (anti-dressing) | NICE | DONE | pigz pipe 2472 vs gzip 3290 variants/s -> parse-bound, GPU not a win |
| S35 | code/s35_frontier.py | AF-dominance frontier (dense M sweep, small n) | MUST | DONE | PC wins only at tiny M (delta +0.008 at M=5k, nsp=500) |
| S36 | code/s36_phased_ibd.py | correct phased IBD (pairwise longest run, self-excluded) | MUST | DONE | planted 0.671, real 0.486 (<=AF); haplotype release channel bounded |
| —   | code/gpu_round_summary.py | merge s27-s36 into one results table | — | DONE | refine-logs/GPU_ROUND_RESULTS.md |

### GPU-round review (run_20261006_350bbd97, pi-review / glm-flash)

| Round | Score | Verdict | Key outcome |
|---|---|---|---|
| R1 | 7.5 | almost | 2 MAJORs: paired "600-seed" was one 200-seed shard (dict.update last-wins); ladder saturation unsupported; + eig-only range slip, undisclosed gauss approx, dangling parser ref, caption/panel mismatches |
| R2 | 8.0 | almost | Both MAJORs fixed + independently re-pooled; eig-only added at biobank scale (0.497–0.511/12 cells); gauss disclosure verified immaterial; residual: Fig 8(a) caption false about its own panel; shard miscount |
| R3 | 8.5 | **ready** | Fig 8(a) now plots refbasis−0.5 (caption true); shard/seed counts match artifact trail; numbers reconciled; merge made first-wins. STOP (converged) |

## Round 11 — limitation closure (2026-10-06) — plan: EXPERIMENT_PLAN.md

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S38 | M0/M1 | per-individual privacy metric on frontier | AF/PC/governed tiers re-expressed | 3-way disjoint | A*, rho, eps_lb, TPR@1%, LR_max | MUST | DONE | pc_loadings A* .610 / af_1dp .187; controls A*≤.04; A*−(2AUC−1)<0 leaking / >0 controls (non-binormal) |
| S39 | M0/M2 | built sole-egress binding | library vs isolated child + audit hook | n/a | bypass matrix, chain verify, governed tier | MUST | DONE | library leak 1.0; isolated natural-egress blocked 1.0; sanctioned-undeclared blocked 1.0; declared flow+chain 1.0; governed af .624 / pc .509 |
| S40 | M0/M2 | joint release-bundle leak | all-objects attacker vs max single | 3-way disjoint | bundle_agg, excess, LOO, structural zeros | MUST | DONE | bundle .958 vs max single .947 (af_site); excess +.011, 5/5 seeds; noise control .839 == af-only |
| S41 | M3 | multi-chrom phased IBD (bp/cM) | chr1/2/22 full-density | per-superpop 3-way | LSR-bp/cM AUROC vs AF; planted ctrl | MUST | DONE | real LSR 0.42–0.48 ≤ AF 0.84–0.87 on all 3 chr; planted 10 cM control 0.70–0.73 |
| S42 | M3 | realistic missingness channel | real dosage + pop-batch + rho | per-superpop 3-way | AUROC vs rho (monotone); rho=0 null | MUST | DONE | rho 0→1: .50→.64→.76→.89→1.00; rho=0 null .50; bonus sweep resolves the magnitude scale |
| S37 | M4 | semi-synthetic biobank cohort | Li-Stephens real-haplotype copying | generated 3-way | s27 keys + LD/relatedness controls | MUST | DONE | real-haplotype cohorts nsp 150→12000: AF leads every cell (delta −0.03…−0.10); in-sample gap →0 with n; LD/relatedness/eig controls match donor |
| S43 | M4 | control reading (caveats) | S37/S41 control battery | — | per-control reading | NICE | DONE | LD-decay, AF-spectrum, eig-fraction match donor at all cells; relatedness within ~10% at 6/7 cells, +11% at nsp150 (141 vs 128, donor std 48 vs cohort 186) |

## Round 12 — residual-limitation closure (2026-10-07) — plan: EXPERIMENT_PLAN_R12.md

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S44 | M0/M1 | per-individual metric on real + cohort | AF/PC channels on 1000G chr22 + s37 cohort | 3-way disjoint | A*, rho, eps_lb(rule-3 floor), TPR@1% | MUST | DONE | real chr22 pc_loadings A* .271 < af_full .350 (ordering preserved); cohort af_full .133 > pc .040; controls A*≤.067; rule-3 eps_lb valid on cohort (δ≥0.01), ordering-only on real |
| S45 | M0/M3 | biobank-corner cohort cell | s37 at nsp=16000, M=3000 | generated 3-way | s27 keys + 4 controls | MUST | DONE | af .529 vs pc .503, delta −.026; refbasis/eigonly .500; controls pass; nsp=33000 needs 27 GB (infeasible on 24 GB) |
| S46 | M0/M2 | second bundle manifest + uncoarsened single | joint attacker vs max single | 3-way disjoint | bundle/max/excess/LOO + noise ctrl | MUST | DONE | manifest 2: bundle .973 vs uncoarsened single .964 (excess +.008, 5/5); coarsened single .947 (+.026); noise control .964==.964 |
| S47 | M0/M0 | native-code egress boundary test | ctypes route vs audit hook | n/a | native_egress_blocked + bypass matrix | MUST | DONE | ctypes.CDLL(None) denied on ctypes.dlopen/dlsym, 3/3; bypass matrix unchanged |

## Round 13 — real-data structure law (2026-10-07) — plan: EXPERIMENT_PLAN_REAL.md

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S48 | M0/M1 | real-data structure ladder | AF vs PC across population cohorts | 3-way disjoint | af_auc, pc_auc, gap, fst_like, refbasis, eigonly | MUST | DONE | `code/s48_real_structure.py` + `s48_scalematched.py`; chr22/1/2 |
| S49 | M1 | real within-pop mechanism + granularity | top-k loadings / eig-only / OOS basis | 3-way disjoint | topk_auc, eig_only, refbasis, oos | MUST | DONE | `code/s49_real_mechanism.py` (0.996–0.999 full, OOS chance) |
| S50 | M2 | multi-chromosome within-pop replication | chr1/2 × GBR/CHB | 3-way disjoint | gap, controls | MUST | DONE | `code/s50_real_multichrom.py` (+0.062…+0.091) |

## Round 14 — codex closure round (2026-10-08/09) — run_id run_20261009_27c9f1a4

Code-level and interpretation closure across three fresh independent Codex threads
(01a11a42 C1 4→C4 7; 01a11b98 C5 5→C7 7; 01a12007 C8 5→C12 7 **ready**). Central
change: the PC-vs-AF comparison now separates attacker post-processing
(`af_conditioned` 0.941) from a curator refit (`pc_conditioned` 0.964). See
`review-stage/AUTO_REVIEW.md` rounds C1–C12 and `CLAIMS_FROM_RESULTS.md` C10–C14.

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S10 | M1 | as-released PC-vs-AF ordering | conditioned AF / curator-refit PC | 3-way disjoint | af_conditioned, pc_conditioned, af_loo_oracle | MUST | DONE | 0.941 / 0.964 / 0.429 (25 seeds); conditioned AF > full-study PC 25/25 |
| S27 | M1 | joint AF+PC estimator refresh | biobank boundary grid | shard | af_only, pc_calibrated, combined | MUST | DONE | shards 0–3 re-run with the joint estimator; gpu_round_summary regenerated |

## Round 15 — limitation closure for the 20-page draft (2026-10-09) — plan: EXPERIMENT_PLAN_LIMITS.md

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| L1 (s51) | M0/M1 | real large-n low-structure ladder | single super-pop pooled (EUR/EAS/AFR) + fixed 3-pop per-group ladder | 3-way disjoint | af_auc, pc_auc, gap, fst_like, controls | MUST | DONE | pooled n=167-220 PC>AF (chr22 +0.044..+0.075, chr1 +0.048..+0.165, chr2 +0.111..+0.181); ladder falls monotonically and crosses on chr22 (+0.120/+0.042/+0.007/−0.026) and chr1 (+0.058/+0.002/−0.010/−0.018); chr2 stays positive, near-parity at the top rung (+0.109/+0.059/+0.018/+0.006) |
| L2 (s38) | M2 | eps rule-of-three floor + inversion example | validate_tradeoff uses the same FNR floor | 3-way disjoint | eps_lb, rule3_valid, inversion | NICE | DONE | rule3 valid=True at every delta on all channels; inversion af_full (AUC .845, eps 4.12) vs af_2dp (.840, 4.19) recorded |
| L3 (s52) | M3 | second haplotype/LD estimator | max-IBS-share + n-close-haps (leave-self-out) | 3-way disjoint | AUROC, planted control | NICE | DONE | real chance (chr22 .488/.496, chr1 .497/.497, chr2 .497/.498); planted 1500-marker copy detected (.798/.797, .776/.788, .772/.779) |
| L4 (s53) | M3 | fixed eigenvalue-informed basis probe | identity / cosine / λ-weighted fixed bases | 3-way disjoint held-out | AUROC | NICE | DONE | all fixed bases at chance (.46-.53) incl λ-weighted on sim + chr22/chr1/chr2 × EUR/EAS/AFR; full loadings .79-.90 |
| L5 (s37) | M2 | fit-fraction sweep at largest cells | n_fit {5k,10k,20k,24k/full} | generated 3-way | gap, pc_auc, af_auc | NICE | DONE | gap stable (5 seeds): nsp8000 −0.034 every n_fit; nsp16000 −0.026 |

## Round 15b — real-gap uncertainty + density axis + seeds (2026-10-11/13)

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S54 | M2 | real-gap permutation null + per-split CI | AF/PC on real s48 cohorts | 3-way disjoint | gap, gap_std, gap_ci, p-value | MUST | DONE | `code/s54_real_permnull.py`; within-pop p<=2e-2, gap CI excludes 0 (GBR +0.085 [+0.058,+0.126], CHB +0.108 [+0.047,+0.202], YRI +0.053 [+0.038,+0.070]); 3pop +0.007 [−0.008,+0.020] near parity/straddles 0; 5pop −0.041 [−0.077,−0.005] |
| S55 | M2 | real within-pop density axis | GBR/CHB/YRI chr22 at 3 densities | 3-way disjoint | af_auc, pc_auc, gap | MUST | DONE | `code/s55_real_density.py`; PC>=AF at ~17k/9.8k/47k markers (gap +0.007/+0.035/+0.000, +0.005/+0.038/+0.000, +0.002/+0.021/+0.000); ordering density-stable |
| S37b | M2 | seed bump 3→5 | semi-synthetic cohort + fitsens | generated 3-way | gap, sd | NICE | DONE | nsp150-16000 gap −0.026…−0.088 (sd up to 0.027); fitsens nsp8000 −0.034 stable, nsp16000 −0.026; merged s37 json rebuilt from the 5-seed shards |
| S42b | M2 | seed bump 3→5 | batch-missingness rho sweep | 3-way disjoint | AUROC vs rho | NICE | DONE | rho 0→1 .490/.629/.755/.893/1.000; candonly .529/.651/.760/.888/1.000 |

## Round 16 — real structure-zero switch + scope closure (2026-10-13) — plan: EXPERIMENT_PLAN_LIMITS.md

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S56 | M0/M1 | real structure-zero switch, fixed-reference design | pooled single-super-pop, study grown + reference held at 150 | 3-way disjoint | af_auc, pc_auc, gap, controls | MUST | DONE | `code/s56_real_poolmax.py`; study grown to n350 with n_ref=150 (AFR full-reference; EUR/EAS capped 147–148): AFR gaps +0.102/+0.140/+0.181 (chr22/1/2, sd 0.018–0.037); zero-structure PC>AF does not flip inside the 1000G pool (EUR 503/EAS 504/AFR 661) — measured bound; per-split gaps exported (superseded the earlier study-heavy run reported in round 15b) |
| S57 | M3 | third haplotype estimator | per-window LD-correlation | 3-way disjoint | AUROC, planted control | NICE | CUT | a within-individual H0/H1 window r² measures phase LD, not membership (no donor/release comparison); s52 already carries two non-LSR statistics — a near-duplicate does not change a reviewer belief |
| S58 | M3 | eigenvalue-dependent attack probe | covariance reconstructed from released eigenvalues | in-simulation | AUROC | NICE | CUT | the paper already states the eigenvalue-only negative is scoped to the tested data-independent bases; a further in-simulation probe is speculative |
| D1 | M2 | ε convention + TCB note | text only | — | — | NICE | TODO | main.tex |

## Round 17 — cross-model defect closure after the round-16 audit (2026-10-14)

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S59 | M0/M1 | GPU AUROC tie-handling disclosure | constant + half-tied + archived cell | n/a | cpu_midrank vs gpu_argsort_rule, tie_mass | MUST | DONE | `code/s59_auroc_tiecheck.py`; constant [1,1,1,1]/[1,1,0,0] → CPU 0.5 vs GPU rule 0.0; half-tied 0.505 vs 0.252; archived-s28-cell tie_mass 0.0 (n=3000) → headline unaffected |
| S60 | M1 | monotone share + per-seed crossing | s51 fixed-composition ladder | 3-way disjoint | cumulative squared-proj share, crossing location | MUST | DONE | `code/s60_share_crossing.py`; share monotone in k on every rung; chr22 crossing mean n≈341 sd≈44 (per-seed: seed0 three crossings, seed4 none; 3 single-crossing seeds) |
| S56 | M1 | s56 fixed-reference rewrite | reference held at 150, non-members shrink | 3-way disjoint | af_auc, pc_auc, gap, controls | MUST | DONE | re-run chr22/1/2: AFR n_ref=150 through n350 (gap +0.102/+0.140/+0.181, sd 0.018–0.037); EUR/EAS n_ref 147–148 at n350 |
| D17a | M1 | s54 scope + gap_p_gt0 label | text only | — | — | MUST | DONE | paper §scale: "five tested cohorts not eighteen"; p tests channel separation not PC−AF |
| D17b | M1 | s37 controls re-statement | text only | — | — | MUST | DONE | AF/eig match all cells; LD matches 5 cells, undefined at M=3000/6000 |
| D17c | M1 | s56 wording corrected | text only | — | — | MUST | DONE | dropped "necessarily falls to 30"; control range 0.37–0.63; n350 AFR quoted |
| D17d | M1 | s42 "identically" corrected | text only | — | — | MUST | DONE | two curves track within ≈0.02 interior rho (candonly marginally higher at rho=0.75) |
| D17e | M2 | CLAIMS stale s56/s54/pending-permutation | text only | — | — | MUST | DONE | C6 unsupported-entry + Next-experiments updated |

## Round 19 — selection-input spectral tomography (deep-search direction, 2026-10-16)

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S61 | M0 | interface witness + preservation audit | real 1000G chr22, AF-preserving block swap | 2 cohorts x 5 seeds | eigen-equation residual, histogram/AF/variance equality, genotype-changed counts | MUST | DONE | `code/s61_spectral_tomography.py`; identity `Y=C_H V` residual 3.8e-11; per-SNP histograms/AF/variance exact; genotypes changed 0-20 per swap; row heterozygosity shifts 1.8e-3-4.9e-3 |
| S61 | M1/M2/M3 | coverage + target + membership vs #outputs | passive/random/adaptive x q{1,2,4,8} | 2 cohorts x 5 seeds | visible_share, scatter_recon_err, target_resid, membership AUROC, perm null | MUST | DONE | balanced conditioned-AF baseline 0.813 (EUR)/0.718 (3pop); q8 random gap vs AFcond +0.186/+0.189 (5/5); adaptive-vs-random CI includes 0 at every q |
| D19a | M2 | experiment plan | text only | — | — | MUST | DONE | `refine-logs/EXPERIMENT_PLAN_DEEP_20261016.md` |
| D19b | M2 | results summary | text only | — | — | MUST | DONE | `refine-logs/EXPERIMENT_RESULTS_DEEP_20261016.md` |

## Round 20 — selection-input spectral tomography (deep review loop, 2026-10-16)

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S62 | M1 | honest-sample-size sensitivity | nested honest n=120 vs 180, partner/ref/shadow fixed | chr22 EUR, 5 splits | gap vs strongest passive | MUST | DONE | `code/s62_sample_size.py`; q8 random +0.015 (n120) vs +0.051 (n180); single release loses at q1 both arms |
| S63 | M1/M2 | real PLINK 2.0 export read-back | `plink2 --pca 10 allele-wts` on real EUR chr22 VCF | chr22 EUR, 5 seeds | identity residual, membership AUC | MUST | DONE | `code/s63_plink_export.py`; identity residual max 6.2e-6 (40 exports); q8 read-back 0.994–1.000; q8 gap vs strongest passive +0.014 [+0.007,+0.021], 5/5 |
| S61 | M0–M3 | selection-input spectral tomography | real 1000G chr22/chr1, AF-preserving swaps | 3 cohorts x 5 seeds | coverage/recovery/membership | MUST | DONE | `code/s61_spectral_tomography.py`; see `refine-logs/EXPERIMENT_RESULTS_DEEP_20261016.md` |
| D20a | M2 | deep review loop (codex) | auto-review-loop, 10 rounds | — | — | MUST | DONE | `review-stage/DEEP_REVIEW_20261016.md`; 4→7/10, verdict ready |

## Round 21 — spectral-tomography round-2 closure (deep loop, pi-review, 2026-10-16)

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S61 | M1 | TPR@5%FPR + clipped no-release control | real 1000G chr22/chr1 | 3 cohorts x 5 seeds | membership_tpr_fpr05, no_release_controls, adaptive_minus_random | MUST | DONE | `code/s61_spectral_tomography.py`; q8 TPR@5% = 1.000/0.580/1.000; single-release 0.675/0.139/0.616; clipped controls 0.455–0.501 |
| S62 | M1 | TPR column added | nested honest n=120/180 | chr22 EUR, 5 splits | membership_tpr_fpr05 | MUST | DONE | `code/s62_sample_size.py`; q8 TPR 1.000/0.995 |
| S63 | M1/M2 | TPR column added | PLINK 2.0 read-back | chr22 EUR, 5 seeds | membership_tpr_fpr05 | MUST | DONE | `code/s63_plink_export.py`; q8 read-back TPR 1.000/0.962 vs baseline 0.833–0.984 |
| S64 | M2 | q-saturation cell | passive/random/adaptive x q{1..16} | 3 cohorts x 5 seeds | gap vs strongest passive at q16 | MUST | DONE | `code/s64_saturation.py`; within-pop saturates at q8 (TPR already 1.000); 3pop keeps gaining +0.049/+0.048 (q8→q16) |
| S65 | M0 | identity ablation control | identity stack vs raw-V stack, q{1,2,4,8} | 3 cohorts x 5 seeds | gap vs strongest passive | MUST | DONE | `code/s65_identity_ablation.py`; raw-V q8 −0.004/−0.012/−0.050 vs identity +0.015/+0.123/+0.013 |
| S66 | M0 | control-provenance check | cal-span AUC under 3 preprocessings | 3 cohorts x 5 seeds | control AUC | MUST | DONE | `code/s66_control_provenance.py`; 27 arm-means in 0.455–0.509 (none exceeds 0.509) across public / own-center-scale / raw-centered |
| D21a | M2 | round-2 integrity re-audit | codex read-only | — | — | MUST | DONE | thread 01a14455; verdict FAIL → WARN → PASS within the declared scope; all reporting defects closed |

## Round 22 — fair comparison and interface confirmation (GPT analysis, 2026-10-17)

| Run ID | Milestone | Purpose | System / Variant | Split | Metrics | Priority | Status | Notes |
|--------|-----------|---------|------------------|-------|---------|----------|--------|-------|
| S67 | M1 | strongest passive adversary bar (gate) | oracle held-out best of af_cond/spec_cond/fused_z/fused_cal | 3 cohorts x 5 seeds x 2 chr | gap_vs_oracle (AUC + TPR@5%FPR), CIs | MUST | DONE | `code/s67_joint_passive.py`; q8 AUC gap +0.0066 EUR / +0.0656 3pop / +0.0037 EAS (chr22), +0.0130/+0.0443/+0.0048 (chr1); TPR gaps larger; 5/5 (chr22 EAS 4/5); center/scale-only 0.69–0.81 |
| S68 | M0/M2 | interface manifest + approx solver (gate) | PLINK 2.0 `--pca approx allele-wts` read-back | chr22 EUR, 5 seeds | identity residual, q8 AUROC, span move | MUST | DONE | `code/s68_interface.py`; manifest all items obtainable except honest matrix; approx = exact eigenvalues/span at n=240; identical residual 6.242e-6; solver-seed span move 6e-8 mean |
| S69 | M1 | precision-reliable directions | reliability floor on the recovered operand Y | 3 cohorts x 5 seeds x 2 chr | geometric vs reliable rank at float64/32/6dp/4dp | MUST | DONE | `code/s69_reliable_dirs.py`; uses Y_t=V Lambda - K V; chr22 EUR 78→30 (6dp), 78→11 (4dp); 3pop reliable AUROC 0.806→0.711 (6dp) |
| S70 | M1 | known-member-fraction sweep | 0/5/10/50% of honest N, drawn from calibration half; mechanism-sign orientation | 3 cohorts x 5 seeds x 2 chr | gap vs label-free fused_z | MUST | DONE | `code/s70_known_fraction.py`; positive 5/5 at every fraction; EUR +0.064 (0%) → +0.011 (50%) |
| S71 | M2 | scientific cost under protection | fixed vs input-varying honest GWAS | 2 cohorts x 5 seeds x 2 chr | lambda_GC, power@1e-4, off-causal fraction, effect agreement, QC | NICE | DONE | `code/s71_protection_cost.py`; one QC marker universe; lambda/power/causal-count/off-causal/agreement unchanged; power low at n=120/360; 3pop off-causal 1.9e-4 |
| S72 | M2 | matched printed-output comparison (closes J9) | q8 PLINK 2.0 `--pca approx allele-wts` read-back vs the passive bar from the same outputs | 2 cohorts (EUR, 3pop) x 5 seeds x 2 chr | gap_q8_vs_oracle_single / _q (AUC + TPR@5%FPR), CIs | MUST | DONE | `code/s72_printed_matched.py`; positive-mean AUC gap every cell (+0.0054 EUR chr22 / +0.0100 EUR chr1, +0.0082 / +0.0060 3pop); interval excludes zero only chr1 EUR; oracle_q ≈ oracle_single |
| D22a | M2 | claims | text only | — | — | MUST | DONE | `CLAIMS_FROM_RESULTS_ARB_20261017.md` |
| D22b | M2 | codex auto-review-loop rounds 1-4 + fixes | codex read-only, thread 01a149b1 | — | — | MUST | DONE | R1 6/10 "No" (matched printed comparison, J7 FP interpretation, J5 scope, per-seed auditability); R2 6/10 "No" (s72 fused_cal trained on conditioned features; all-q interpretation; cohort path collision; TPR count); R3 7/10 "Almost" (headline qualification, TPR counts, J7 wording, two-cohort tracker); R4 7/10 **READY**. All fixes applied: s72 fused_cal retrained on unconditioned features, s72 info model + late-fusion labeling, cohort-tagged paths, per_split + ci_definition, metric-specific TPR oracle, split-count/interval wording |


