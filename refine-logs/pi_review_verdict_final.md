# Final submission verdict — round-5 fix verification (2026-10-03)

Verification scope: every number in `paper/main.tex` re-checked against the
fresh artifacts (`results/raw/*.json`, all post-dating the maf fix); the s7/s8
script fixes verified in code; fedlmm2025's corrected authors externally
verified (Xiaowen Suo, Fuzhong Xue, Yanyan Zhao — bioRxiv 2025.12.16.693409).

## Fix verification — all four applied

1. **Unified panels (n_snps=4000)** — YES. s1c/s5 rebuilt on 4000 SNPs; Table 1
   matches s1c exactly (0.590/0.622/0.840/0.845/0.964, TPRs and CIs, control
   0.502, per-ancestor 0.860/0.849/0.851); boundary table matches s5
   (0.270/0.993 … 0.015/0.038/0.042); Table 2 matches s2b (3 seeds: 0.844 /
   0.942 / +0.098) with an honest caption referencing the 5-seed confirmation
   (s2c base_M2 0.9408); the round-4 missingness rounding nit is resolved.
2. **s8 reframed** — YES. Paper section "Robustness across functional forms";
   limitations state "not three independent attack channels… functional-form
   robustness rather than cross-channel generality; a genuinely distinct
   channel (haplotype-based) is future work" — exactly the required downgrade.
   `s8_families.json` now has a config block (seeds, n_site, n_per_group,
   n_snps, control, maf_filter [0.05,0.95]).
3. **S7 fixed and run** — YES. Script verified: full-split GT parsing with both
   alleles ({0,1,2}; the docstring documents the old bug), shuffled stratified
   3-way split, 15 kb distance-thinning. `s7_real.json` matches the paper's
   real-data paragraph exactly: 0.699 full / 0.696 2dp / 0.577 1dp / 0.646
   top-200 / 0.793 per-population; reference-basis structure control 0.510
   (chance). Meta: 1363 SNPs, 450/panel.
4. **Paper regenerated** — YES. The mechanism paragraph now correctly explains
   the sign (members farther in simulation; naive orientation 0.12 / flipped
   0.88, matching s1c study_pc_distance 0.118 and s2c pcdist 0.881); the
   ancestry control paragraph describes the corrected reference-fitted basis;
   the dimension sweep now matches s6 including the saturation corner; the
   bibliography fixes landed (Berrang et al.(2026); Pascoal et al.(2023);
   Suo/Xue/Zhao, externally verified).

Every table and in-text figure traces to a fresh artifact — with one exception
(§ blocking item below).

## Answers

**(1) Internally consistent?** Almost — one exception. All numbers are fresh
and artifact-backed except the real-data paragraph's claim that the
PC-diagnostic leak "was not separately evaluable on this real release":
`s7_real.json` contains `study_pc_distance_auroc = 0.343` (naive orientation).
The number exists; it is below chance — on real data members sit *closer* to
the study centroid than non-members, reversing the simulated sign (consistent
with the leverage effect dominating direction-tuning at this operating point;
flipped orientation ≈ 0.657). "Not separately evaluable" contradicts the
artifact, and a reviewer with the artifacts reads it as selective omission.

**(2) Submission-plausible as a measurement paper?** YES. The paper now sits on
a single coherent pipeline (fixed maf, unified panels), every table matches its
artifact, the mechanism is properly signed and controlled (independent-spectrum
collapse +0.003, attacker-realistic centering intact +0.10±0.01, permutation
0.486, eigenvalue ratio 1.000), the real-data paragraph is honest about
magnitude (ordering reproduces, smaller), the functional-form reframing is
exactly right, and the corrected pipeline *strengthened* the result (leakage
0.84–0.96 rather than 0.72–0.96; the real-data ordering reproduces at 0.58–0.79).

**(3) Remaining blocking item:** one, text-only — **report the real-data
study-PC-distance number honestly.** Replace "not separately evaluable" with
the measured fact: 0.343 naive / ≈0.657 flipped, sign reversed relative to the
simulation, weak and not interpreted as evidence for or against the PC channel
on real data. Everything else is polish:

- T3 range: paper says "+0.10 to +0.12"; per-seed is 0.092–0.116 — say
  "+0.09 to +0.12, mean +0.10".
- The abstract's "inflates ≈0.2" derives from the old-pipeline biased-vs-
  disjoint comparison (0.916 vs 0.720). Either re-run the biased arm with the
  fixed pipeline (`split_population3` with n_nonmem=0 reproduces it — cheap) or
  footnote that the inflation figure predates the maf fix.
- prov-o bibitem still reads "Lebeck & Gill(2022)" for the W3C PROV-O
  Recommendation (2013; editors Moreau & Missier).
- State the 15 kb thinning window in the real-data paragraph.

## One-line PC summary

A controlled measurement paper showing the released study-fitted PC spectrum —
not the summary statistic — is the dominant intermediate membership-leakage
channel in multi-site GWAS (0.88 alone, 0.94 combined, vs 0.84–0.85 for AF
releases; per-site release 0.96), with the mechanism identified as study-fit
in-sample-ness under clean disjoint controls, functional-form robustness,
an honest real-data check on 1000 Genomes chr22 reproducing the ordering at
smaller magnitude, and plainly stated limits — submission-ready once the
real-data PC-distance number is reported rather than omitted.
