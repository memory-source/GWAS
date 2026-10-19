# Cross-model review — real-data design + second attack family

You are the cross-family reviewer. Two round-3/4 critiques are being addressed
in this round; check the approach before it is finalized.

## (A) Second/third attack family ('one working family' critique)
Replaced the structurally-weak beacon attack with two structurally distinct
families, run on the SAME panels as the AF-LRT (5 seeds):
- F2 AF-correlation: per-individual corr(genotype, released study-AF vector) −
  corr(genotype, reference-AF vector).
- F3 ML-deviation: cross-validated logistic regression on per-individual mean
  squared deviation from study-AF vs reference-AF.
Fresh results (results/raw/s8_families.json):
  F5: AF-LRT 0.845 / AF-corr 0.833 / ML-dev 0.832 / PC-distance 0.882
  F3(2dp): 0.840 / 0.830 / 0.829 / 0.882
  F4(1dp): 0.622 / 0.660 / 0.659 / 0.882
So three structurally different families agree on the frontier ordering, and the
PC-diagnostic channel still overtops all of them.

## (B) Real-data external validity (S7)
Plan: download the public 1000 Genomes Phase-3 chr22 VCF + integrated_call
panel; take three super-population pseudo-sites (EUR/EAS/AFR); split each three
ways disjointly (study / attacker reference / independent non-members); MAF-filter
[0.05,0.95]; LD-thin (~100 kb window); rerun the frontier tiers and the
PC-spectrum mechanism. No simulator assumption in the control (ancestry-matched
by construction). Script: code/s7_real.py.

## Questions
1. Are F2/F3 legitimate distinct membership-inference families, or are they
   re-parameterizations of F1? If the latter, say so plainly.
2. Any correctness risk in the S7 design (leakage of the split, LD, MAF,
   per-population AF release)?
3. Does resolving (A)+(B) close the 'all simulation / one attack family'
   rejection clauses, or is something still missing? Be concrete; do not invent
   requirements.
Read /data/lvvyongji1/GWAS/code/s7_real.py and results/raw/s8_families.json.
