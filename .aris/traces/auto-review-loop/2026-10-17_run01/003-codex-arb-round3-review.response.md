**Score: 7/10. Verdict: Almost.**

The substantive round-2 defects are fixed. I found **no remaining critical numerical or implementation defect** in the revised s72 comparison. The narrowed result is now honestly presented: float64 superiority is established against the four specified constructions; printed-output evidence is limited, with a positive conditional AUC interval only for chr1 EUR.

Verification supports that assessment:

- All 16 saved s72 AUC/TPR gap summaries—including bootstrap endpoints and positive counts—reconstruct exactly from `per_split`.
- All **160 expected cohort-specific exports** are preserved, with completed PLINK logs.
- A direct chr1 EUR seed0 replay from eight printed exports exactly reproduced active AUROC **0.9968443** and conditioned-single AUROC **0.9828191**.
- Unconditioned calibration fitting is restored correctly ([s72:141–154](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:141)); the information model and late-fusion qualification are now accurate ([s72:114–119](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:114)).

### Remaining weaknesses and minimum fixes

1. **The practical result remains limited; the headline needs one final qualification.**

   The corrected chr1 EUR gap is **+0.01005 [0.00380, 0.01645]** ([raw chr1:151–164](/data/lvvyongji1/GWAS/results/raw/s72_printed_matched_chr1.json:151)). Chr22 EUR and both 3pop AUC intervals cross zero. That is correctly disclosed.

   However, the center claim ends with “the interpretable effect is cross-population” ([claims:15–16](/data/lvvyongji1/GWAS/CLAIMS_FROM_RESULTS_ARB_20261017.md:15)). The larger supported cross-population increment belongs to the **float64** experiment; printed cross-population superiority remains unresolved ([raw chr22:318–337](/data/lvvyongji1/GWAS/results/raw/s72_printed_matched_chr22.json:318), [raw chr1:318–337](/data/lvvyongji1/GWAS/results/raw/s72_printed_matched_chr1.json:318)).

   **Minimum fix:** write: “The larger float64 increment is cross-population; matched printed-output evidence is positive with an interval excluding zero only for chr1 EUR.” No additional experiment is required if this remains the submission’s scope.

2. **The results document still contradicts the corrected TPR counts and interval.**

   [Results:59–60](/data/lvvyongji1/GWAS/refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md:59) incorrectly says all chr22 TPR splits are positive and that chr1 EAS’s TPR interval crosses zero. [Results:264](/data/lvvyongji1/GWAS/refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md:264) still says “positive 5/5 within populations.”

   The actual EAS counts are **4/5 on chr22** and **3/5 on chr1**. Chr1 EAS’s interval is **[+0.003333, +0.029796]**, which excludes zero ([raw chr22:2030–2035](/data/lvvyongji1/GWAS/results/raw/s67_joint_passive_chr22.json:2030), [raw chr1:2029–2034](/data/lvvyongji1/GWAS/results/raw/s67_joint_passive_chr1.json:2029)). The claims file now gets these counts right.

   **Minimum fix:** reconcile those three sentences with the raw values. These are reporting errors, not failures of the calculations.

3. **J7 retains the unsupported nominal-error and causal-attribution wording in two places.**

   “Above nominal, partly LD proxies” remains in [claims:38](/data/lvvyongji1/GWAS/CLAIMS_FROM_RESULTS_ARB_20261017.md:38), and the same interpretation survives in [results:303–305](/data/lvvyongji1/GWAS/refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md:303).

   The implemented quantity counts significant unplanted markers ([s71:111–125](/data/lvvyongji1/GWAS/code/s71_protection_cost.py:111)); its saved value is correctly computed ([raw:57–59](/data/lvvyongji1/GWAS/results/raw/s71_protection_cost_chr22.json:57)). It has no clean nominal Type-I-error reference, and the experiment does not quantify how much LD caused.

   **Minimum fix:** replace both passages with: “The off-causal discovery fraction is \(1.9\times10^{-4}\); causal proxies can contribute, so it is not a null-error estimate.”

The revised QC statement is correct when “per-marker scatter” means the per-marker sum of squares. Columnwise genotype permutations preserve that quantity; cross-marker scatter changes, enabling the mechanism. Also correct the tracker’s s72 count from **three cohorts to two** ([tracker:220](/data/lvvyongji1/GWAS/refine-logs/EXPERIMENT_TRACKER.md:220)); an additional EAS experiment is unnecessary.

**The remaining required changes are text-only.** I would not require another large experiment for this bounded mechanism-and-measurement submission.

### Memory Update

- **Closed:** calibrated-comparator regression, access-model inconsistency, export collisions and numerical auditability.
- **Confirmed:** J9 is partial; only chr1 EUR has a positive conditional AUC interval.
- **Unresolved boundary:** broader printed-output and independent-subject generalization.
- **Pattern:** implementation repairs are complete, but stale claims persist in secondary summary paragraphs.