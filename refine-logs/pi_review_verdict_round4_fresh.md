# Senior Review — Round 4 (fresh run, targets: Limitations / replication)

Reviewer: independent artifact-level audit (paper + code + raw JSONs read in full; ~45 paper
numbers matched against `results/raw/*.json`; caches and compile checked). Prior rounds:
7.0 → 7.5 → 8.0 "ready" (run_20261004_bcd63f1a). This round reviews the new §6
(Replication and a Second Channel), Fig. 6, and the rewritten Limitations.

## What holds up (verified, not trusted)

- **Numbers trace 1:1 to artifacts.** All s16 values (chr1 5257 / chr2 4664 markers; AF
  0.768±0.045 / 0.752±0.046; PC 0.675 / 0.689; AF+PC 0.805 / 0.799; mechanism in ~0.73–0.74
  vs oos/swapped 0.51; granularity full 0.675/0.690, eig-only 0.516/0.514), all s18 values
  (PC 0.691±0.013 vs AF 0.627±0.015; top1 0.524 / top5 0.621 / eig 0.512; A_in 0.967,
  A_oos 0.547), s19 (hap 0.990, OOS 0.487, AF 0.845), s17 (n1200: 0.909/0.896/0.966;
  n3000: 0.807/0.745/0.862), and s20 (M3 sweep 0.49/0.69/0.93/1.00; batch gloss 0.76 backed)
  match the JSONs exactly.
- **s16 mechanism design is clean**: A/B/C are same-ancestry-composition disjoint splits
  (permute-within-superpopulation), so own/other/swapped are valid in-ancestry nulls.
- **s17 consumes the released object** (`rel["loadings"]` with mask), train-half orientation
  calibration, out-of-fold logistic combination — consistent with the W1 fix.
- **s18 is a genuine recombination-realistic cross-check** (real recombination 1.25e-8,
  5 Mb, realistic SFS), not the injected-LD parametric model.
- **s19's OOS transfer control is sound**: both pair tables fit on panels that saw neither
  study nor targets; the in-sample condition is the honest self-fit.
- **s20's M3 channel is a genuine missingness channel**: released per-SNP call-rate table
  estimated out-of-sample (test members excluded), null control at equal rates (0.506),
  self-ref variant, graded sweep. The Discussion's 0.76 batch gloss now cites the
  artifact-backed s15 numbers (0.846 / 0.800 / 0.759) — prior error fixed.
- **Real data is real**: chr1/chr2 1000G Phase-3 VCFs (1.1/1.2 GB) present with correct
  release naming; caches match the paper's marker counts (5257/4664; 1668 EUR/EAS/AFR
  samples); MAF [0.05,0.95] + 15 kb thinning applied in the loader.
- **New honest boundary reported**: PC-dominance does not survive at n≥1200 (AF 0.909 vs
  PC 0.896) — disclosed, not tuned away. Limitations rewritten accordingly.
- Compiles fresh: 10 pages, 0 errors, 0 undefined refs (3 cosmetic overfull boxes). Fig. 6
  is generated from s16/s17/s18/s19 JSONs.

## 1. Score: 8/10

## 2. Critical weaknesses (ranked)

1. **Internal ordering-claim contradiction (abstract ×2, contribution 1, Discussion) — the
   previously-fixed error class regressed under the new content.** The abstract claims
   "the PC diagnostic is the largest intermediate-release leakage channel on the simulated
   panels" and "per-site AF release is the largest single channel of all (0.96)", yet the
   same abstract/paper reports the hap-LD channel at **0.99** and the membership-correlated
   M3 missingness channel at up to **1.00** — both on simulated panels, both larger than
   PC (0.88) and per-site AF (0.96). The intended scoping (canonical workflow diagnostics;
   marginal-frequency families) is never stated.
2. **s18's "swapped" control is structure-contaminated and reported without dispersion.**
   In s18 the "sites" are different msprime populations (A, B, C = populations 0/1/2), so
   A_swapped (population B vs C under A's basis) is not a guaranteed-null: per-seed
   0.704/0.516/0.596/0.504/0.509 → mean 0.566 ± 0.077, and seed 0 is 0.70. The paper
   reports "swapped 0.566" as if a clean at-chance null, and calls the design "the same
   held-out-site mechanism" although (unlike s14/s16) the non-member set is a different
   population, so A_in 0.967 also includes ancestry separation. The mechanism claim still
   stands on A_in 0.967 vs A_oos 0.547 (ancestry separation under an out-of-sample basis is
   negligible) — but that argument is not made.
3. **"+0.056 at both" (§6 site-size paragraph)**: s17 gives +0.0565 (n1200) and +0.0550
   (n3000); the second rounds to 0.055.
4. **s17 docstring overpromises vs what ran** (claims n∈{1200,…,12000}, M=20000, 20 seeds,
   LD-dose range; ran {1200,3000}, M=10000, 8 seeds, ld100+ld0 — the paper matches the
   run, the docstring does not); §4's seed inventory does not cover s16–s20 (s18: 5 seeds;
   s19: 5; s20: 10/5/3).
5. **s20's structure block is uncited and unreconciled**: chr22 Fst-like 0.0725 (first-2000-
   marker subset) vs the paper's s15-based ≈0.10 (full chromosome) — a cross-script
   discrepancy a careful reader will find; also chr1 0.142 / chr2 0.114 are *more*
   structured than chr22, which would actually strengthen the shared-structure hypothesis
   if cited.
6. **Hap-LD operating point and release realism unstated**: 900/900/900, 4000 SNPs, 5
   seeds, parametric sim; the "released joint LD table" is a hypothetical release object —
   one clause (imputation/LD reference panels do ship fitted LD) would tighten it.

## 3. Minimum fixes

- **W1 (must-fix, text only)**: re-scope the three "largest/dominant" claims — abstract:
  "the largest of the canonical GWAS diagnostic releases (λGC, PC spectrum, i.i.d.-model
  missingness, condition bit) on the simulated panels"; "per-site AF ... the largest
  marginal-frequency single channel (0.96)"; contribution 1 heading likewise; Discussion
  sentence likewise — and state explicitly that the hap-LD (0.99) and M3-batch (≤1.00)
  channels exceed PC.
- **W2**: report "swapped 0.566 ± 0.077 (seed 0: 0.70)" plus one sentence disclosing that
  the msprime mechanism uses populations as sites (so the swapped statistic is
  structure-sensitive) and that the clean at-chance control is the out-of-sample basis
  (0.547 ± 0.020). Better (optional, ~10 lines + no new msprime sim): rerun the s18
  mechanism with within-population three-way splits.
- **W3**: "+0.056/+0.055" or "+0.055–0.056".
- **W4**: fix the s17 docstring; extend the §4 seed inventory with s16–s20 counts.
- **W5**: one sentence reconciling the two chr22 structure estimates (marker subset) and
  citing chr1/chr2 as more structured.
- **W6**: one sentence stating the hap-LD operating point and the release-realism clause.

## 4. READY? Almost

All measurements verify; the remaining items are text-level (~an hour) plus one optional
cheap rerun. W1 is a must-fix before submission: the headline ordering claim is
contradicted by the paper's own new numbers in the same abstract.
