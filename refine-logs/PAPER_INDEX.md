# Paper index — auditable reveal governance for multi-site GWAS

Maps every paper artifact to its generating code and raw results. Canonical
numbers live in `results/raw/*.json`; paper text is in `paper/main.tex`.

## Main-text sections → experiments → results

| Paper section | Experiments | Results | Figures |
|---|---|---|---|
| Controlled MI protocol | s0_sanity, s0b_sanity_seeds, common.py | `results/raw/s0*.json` | — |
| The diagnostics release dominates | s1c_frontier, s1d_bias, s2b_intermediate, s2c_mechanism | `s1c_frontier.json`, `s1d_bias.json`, `s2b_intermediate.json`, `s2c_mechanism.json` | fig1, fig2, fig5 |
| Site-size and precision boundary | s3b_paradox, s4_msprime, s5_boundary, s6_dim | `s3b_paradox.json`, `s4_msprime.json`, `s5_boundary.json`, `s6_dim*.json` | fig3, fig4 |
| Replication and a second channel | s7_real, s8_two_family, s9_orientation, s10_paired, s11_robust, s16_multichrom, s17_scale, s18_msprime_pc, s19_haplotype | `s7_real.json`, `s16_multichrom.json`, `s17_scale.json`, `s18_msprime_pc.json`, `s19_haplotype.json` | fig6 |
| Building and stress-testing the governance layer | governance.py, s21_sensitivity, s22_admixture, s23_governed, s24_crossrun, s25_paradox, real_phased.py, s26_ibd | `s21_sensitivity.json`, `s22_admixture.json`, `s23_governed.json`, `s24_crossrun.json`, `s25_paradox.json`, `s26_ibd.json` | fig7 |
| **Scale, density, ancestry (GPU round)** | common_gpu.py, s27_gpu_scale, s28_permnull, s29_negation, s29_ladder, s30_density, s31_fivepop, s33_meta_attack, s34_parse_bench, s35_frontier, s36_phased_ibd | `s27_gpu_scale_shard*.json`, `s28_permnull_shard*.json`, `s29_negation_shard*.json`, `s29_ladder.json`, `s30_density.json`, `s31_fivepop.json`, `s33_meta_attack.json`, `s34_parse_bench.json`, `s35_frontier.json`, `s36_phased_ibd.json` | fig8 |
| Discussion / limitations | (aggregates the above) | `gpu_round_summary.json` | — |

## Result files → code (1:1)

- `s1c_frontier.json` ← `code/s1c_frontier.py`
- `s1d_bias.json` ← `code/s1d_bias.py`
- `s2b_intermediate.json` ← `code/s2b_intermediate.py`
- `s2c_mechanism.json` ← `code/s2c_mechanism.py`
- `s3b_paradox.json` ← `code/s3b_paradox.py`
- `s4_msprime.json` ← `code/s4_msprime.py`
- `s5_boundary.json` ← `code/s5_boundary.py`
- `s6_dim*.json` ← `code/s6_dim.py`
- `s7_real.json` ← `code/s7_real.py`
- `s8_two_family.json` ← `code/s8_two_family.py`
- `s9_orientation.json` ← `code/s9_orientation.py`
- `s10_paired.json` ← `code/s10_paired.py`
- `s11_robust.json` ← `code/s11_robust.py`
- `s12_crossover.json` ← `code/s12_crossover.py`
- `s13_real_loadings.json` ← `code/s13_real_loadings.py`
- `s14_real_mechanism.json` ← `code/s14_real_mechanism.py`
- `s15_controls.json` ← `code/s15_controls.py`
- `s16_multichrom.json` ← `code/s16_multichrom.py`
- `s17_scale.json` ← `code/s17_scale.py`
- `s18_msprime_pc.json` ← `code/s18_msprime_pc.py`
- `s19_haplotype.json` ← `code/s19_haplotype.py`
- `s20_extra.json` ← `code/s20_extra.py`
- `s21_sensitivity.json` ← `code/s21_sensitivity.py`
- `s22_admixture.json` ← `code/s22_admixture.py`
- `s23_governed.json` ← `code/s23_governed.py`
- `s24_crossrun.json` ← `code/s24_crossrun.py`
- `s25_paradox.json` ← `code/s25_paradox.py`
- `s26_ibd.json` ← `code/s26_ibd.py` (+ `code/real_phased.py`)
- `s27_gpu_scale_shard{0,1,2}.json` ← `code/s27_gpu_scale.py` (+ `code/common_gpu.py`)
- `s28_permnull_shard{0,1}.json` ← `code/s28_permnull.py`
- `s29_negation_shard{0,1,2}.json` ← `code/s29_negation.py` (paired block; pooled by `code/gpu_round_summary.py`)
- `s29_ladder.json` ← `code/s29_ladder.py`
- `s30_density.json` ← `code/s30_density.py` (+ `code/real_multi.py`)
- `s31_fivepop.json` ← `code/s31_fivepop.py` (+ `code/real_multi.py`)
- `s32_phased_ibd.json` ← `code/s32_phased_ibd.py` (superseded by s36)
- `s33_meta_attack.json` ← `code/s33_meta_attack.py`
- `s34_parse_bench.json` ← `code/s34_parse_bench.py`
- `s35_frontier.json` ← `code/s35_frontier.py`
- `s36_phased_ibd.json` ← `code/s36_phased_ibd.py`
- `gpu_round_summary.json` ← `code/gpu_round_summary.py` (consolidated GPU-round table)

## Shared data modules

- `code/common.py` — CPU simulator, `release_objects`, AF-LRT/PC/beacon attacks, AUROC, PCA.
- `code/common_gpu.py` — torch port: host-resident simulator (beta & gauss AF draws), streamed
  AF-LRT, projection distance, blocked exact-Gram PCA (TF32 off), rank AUROC.
- `code/real_panel.py` — cached 1000G chr22 parse.
- `code/real_multi.py` — chr1/chr2 loader, 3-pop and 5-pop.
- `code/real_phased.py` — phased H0/H1 loader (keeps 1000G phase).
- `code/governance.py` — reveal manifest, egress monitor, hash-chained ledger.

## Data

- `data/external/chr22.vcf.gz` — 1000G Phase-3 chr22.
- `data/external/ALL.chr{1,2}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`.
- `data/external/1000g_panel.txt` — sample → population → super-population.
- `results/cache/*.npz` — parsed panels (chr22/chr1/chr2, 3-pop and 5-pop variants).

## Environments

- GPU / torch kernels: `/opt/miniforge/bin/python3` (torch 2.2.0+cu121).
- CPU analyses: `/data/lvvyongji1/GWAS/env/bin/python` (numpy 1.26, scipy, sklearn, msprime).
