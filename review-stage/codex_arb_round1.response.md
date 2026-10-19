**Score: 6/10. READY for submission: No.**

The narrowing is honest and makes **J1 evidence-aligned**: the float64 attack exceeds the held-out AUROC oracle over the four specified passive constructions. It does not establish superiority over arbitrary passive processing, and the current claim correctly avoids that assertion.

I independently reconstructed all 30 s67 random-q8 comparisons; the aggregate AUROCs, paired gaps and bootstrap endpoints matched the saved JSONs exactly. The full chr22 s71 replay matched all 53 saved numeric fields within \(3.2\times10^{-15}\). I found no calculation error invalidating J1.

The remaining weaknesses, ranked by submission impact, are:

1. **The practical attack result still lacks the decisive matched comparison.**

   s67 measures the four-construction comparison using float64 eigenpairs; s68 computes printed-output exact/approximate attack performance without constructing that passive comparator. See [s67_joint_passive.py:149–205](/data/lvvyongji1/GWAS/code/s67_joint_passive.py:149), [s68_interface.py:154–179](/data/lvvyongji1/GWAS/code/s68_interface.py:154), and the persisted [s68 results:104–109](/data/lvvyongji1/GWAS/results/raw/s68_interface_chr22.json:104).

   The declaration that this conjunction is “NOT established” is correct ([claims:6–15](/data/lvvyongji1/GWAS/CLAIMS_FROM_RESULTS_ARB_20261017.md:6)). This is an evidence gap, not misrepresentation. Nevertheless, it leaves the strongest incremental attack result detached from the actual export interface—particularly in 3pop, where the effect has useful headroom.

   **Minimum fix:** run the existing chr22 3pop splits through matched q1/q8 PLINK exports and compare the accumulated attack against all four passive constructions formed from the same printed outputs. Report the outcome, including failure. A biobank-scale experiment is unnecessary for this immediate question.

2. **J7’s “false-positive” interpretation is unsupported, and “usable” overstates the utility result.**

   The phenotype depends on eight real-genotype markers, but every tested unplanted marker enters the false-positive denominator ([s71_protection_cost.py:62–72,108–122](/data/lvvyongji1/GWAS/code/s71_protection_cost.py:62)). Unplanted markers correlated with causal markers are not clean marginal nulls.

   This concern occurs in the actual experiment. A targeted replay found all four chr22 off-causal discoveries correlated with planted causal markers after the same PC adjustment; maximum absolute correlations ranged from 0.449 to 0.750. These correlations do not prove each discovery’s cause, but they prevent interpreting the count as clean null calibration. The saved 3pop rate is correctly calculated as approximately \(1.925\times10^{-4}\) ([raw results:57–59](/data/lvvyongji1/GWAS/results/raw/s71_protection_cost_chr22.json:57)); attributing its excess to residual population structure is unsupported ([results summary:185–187](/data/lvvyongji1/GWAS/refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md:185)).

   Likewise, power of 0.025/0.050—recovering 0.2/0.4 of eight causal markers—supports low-power summary stability, not demonstrated useful power ([raw results:15–18,53–56](/data/lvvyongji1/GWAS/results/raw/s71_protection_cost_chr22.json:15)).

   **Minimum fix:** rename `fp_rate` to “off-causal discovery fraction,” remove the nominal-calibration explanation, and replace “usable” in J7 with “association-summary stability in this low-power simulation.” Add null/permuted phenotypes only if retaining an actual Type-I-error claim.

3. **J5 remains stronger than its precision experiment warrants.**

   s69 applies `np.round(A, 6)` or `np.round(A, 4)` to NumPy eigenpairs ([s69_reliable_dirs.py:51–59](/data/lvvyongji1/GWAS/code/s69_reliable_dirs.py:51)). PLINK prints six **significant digits** of differently scaled allele weights ([interface manifest:26–29](/data/lvvyongji1/GWAS/results/raw/s68_interface_chr22.json:26)). These quantization experiments are not interchangeable.

   Furthermore, the floor uses the realized difference from an unrounded float64 reference ([s69_reliable_dirs.py:97–110](/data/lvvyongji1/GWAS/code/s69_reliable_dirs.py:97)). It is an offline heuristic. The reported EUR rank reduction 78→30 and 3pop AUROC reduction 0.806→0.711 are valid outputs of that heuristic ([raw results:30–40,80–90](/data/lvvyongji1/GWAS/results/raw/s69_reliable_dirs_chr22.json:30)); they do not establish counts of reliably usable directions from actual PLINK files.

   The explicit “not a certified usable-rank theorem” caveat is correct. The affirmative wording “bounded by the export precision” and “removes spurious directions” still exceeds what was demonstrated.

   **Minimum fix:** describe J5 as **heuristic retained rank under simulated quantization**, with an offline float64 reference. Use actual parsed export operands if retaining a PLINK-specific precision-boundary claim.

4. **The files called “raw results” discard the measurements needed to audit the intervals.**

   s67 creates per-seed `rows`, then persists only aggregate `blk`; s70 does the same ([s67_joint_passive.py:219–249](/data/lvvyongji1/GWAS/code/s67_joint_passive.py:219), [s70_known_fraction.py:173–191](/data/lvvyongji1/GWAS/code/s70_known_fraction.py:173)). Consequently, the saved positive intervals and split counts—for example, [s67 chr22:355–369](/data/lvvyongji1/GWAS/results/raw/s67_joint_passive_chr22.json:355)—cannot be reconstructed from those JSONs alone.

   My reconstruction confirms the central s67 numbers, so this is an auditability weakness rather than evidence of erroneous intervals. Those intervals describe variation over five splits of this fixed panel; they do not establish independent-subject generalization.

   **Minimum fix:** persist the existing per-seed candidate metrics, selected comparator, active metrics and paired differences. Explicitly label the intervals as conditional split uncertainty. No new experiment is required.

Several things are correct and should remain credited. The genotype exchange preserves per-SNP counts; role partitions are disjoint within each split; attack orientation uses calibration labels, not evaluation labels. The randomized PLINK path genuinely executes, and its agreement with the exact solver at this size is accurately reported. J6 correctly uses a fixed mechanism sign and explicitly compares against the label-free fusion.

The TPR calculation is also correct **for the AUROC-selected comparator** ([s67_joint_passive.py:152–158](/data/lvvyongji1/GWAS/code/s67_joint_passive.py:152)). My additional metric-specific oracle calculation still gives a positive chr22 3pop TPR increment: **+0.1436, CI [0.0771, 0.2100]**. Chr1 3pop remains inconclusive. This comparator distinction needs precise labeling, but it does not overturn J1.

## Memory Update

- **Confirmed:** the narrowed float64 four-construction AUROC result reproduces; no central arithmetic failure was found.
- **Unresolved:** superiority over the matched passive comparator through actual printed outputs, especially in 3pop.
- **Interpretation risks:** off-causal discoveries are not necessarily false positives; low-power stability does not demonstrate useful power; heuristic rank reduction does not certify recoverable directions.
- **Review boundary:** the earlier integrity PASS and readiness verdict concern their stated comparisons. They do not resolve the stronger, still-unmeasured conjunction.