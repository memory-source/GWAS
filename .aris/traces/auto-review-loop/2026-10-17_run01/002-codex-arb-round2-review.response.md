**Score: 6/10. Verdict: No.**

There are genuine improvements: J5 is now adequately scoped, per-split persistence closes the auditability concern, and the metric-specific TPR oracle is implemented correctly. I reconstructed 240 s67/s70 gap summaries from the saved rows; all means, bootstrap endpoints and positive counts matched exactly.

**W1 remains open because s72 changes the calibrated passive baseline. Restoring the original construction removes the positive chr22 EUR confidence interval.** The float64 J1 result remains supported.

### Remaining weaknesses and minimum fixes

1. **Critical: s72 weakens the calibrated comparator rather than transporting s67’s four constructions to printed output.**

   s72 fits logistic weights using conditioned calibration features: `fit(Zc[h], y[h])` ([s72_printed_matched.py:143–148](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:143)). Those known members have already been subtracted. In contrast, s67 deliberately trains on **unconditioned** calibration scores and applies those weights to conditioned scores ([s67_joint_passive.py:124–147](/data/lvvyongji1/GWAS/code/s67_joint_passive.py:124)).

   The resulting failure is visible in the raw JSON: calibrated-fusion AUROC is **0.13796** for EUR seed3 and **0.18303** for 3pop seed2 ([raw chr22:110–114](/data/lvvyongji1/GWAS/results/raw/s72_printed_matched_chr22.json:110), [247–251](/data/lvvyongji1/GWAS/results/raw/s72_printed_matched_chr22.json:247)).

   I restored s67’s recipe using retained EUR PLINK exports for identical inputs and solver seeds. All five active AUROCs and the three nonlearned candidate AUROCs matched s72 exactly. The corrected chr22 EUR gap is:

   **+0.005429; conditional bootstrap CI [−0.000277, +0.011013]; 4/5 positive.**

   Thus the advertised positive interval does not survive this correction. Corrected chr22 3pop also remains inconclusive: **+0.008229 [−0.008348, +0.028693]**.

   **Minimum fix:** restore the original calibration procedure, recompute both chromosomes’ passive bars and intervals, and regenerate J9’s claims. Reuse existing exports wherever available. Do not retain the current “positive CI on EUR” statement without that recomputation.

2. **Major: the all-q interpretation does not follow from the executed information model.**

   The code says the passive adversary does **not** know partner scatter ([s72:114–118](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:114)). Its passive functions nevertheless receive `Ys`, already constructed by subtracting that scatter ([s72:101–103](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:101), [154–157](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:154)). The numerical comparison therefore uses partner information, despite the outsider narrative.

   Furthermore, `oracle_q` sums standardized **scalar per-release scores** ([s72:139–148](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:139)). It does not test general processing of the transcript. Under the information actually supplied, concatenating those same `Ys` and invoking the existing active scorer reproduces the active score ([s72:161–163](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:161)).

   The near-equality of the saved bars is real ([raw chr22:315–317](/data/lvvyongji1/GWAS/results/raw/s72_printed_matched_chr22.json:315)); the broader conclusion about handing over all q outputs is unsupported.

   **Minimum fix:** retain the original same-partner information model, correct the access description, and label `oracle_q` as a diagnostic over the specified **late-fusion constructions**. Remove the general “all q outputs do not strengthen the passive adversary” conclusion. No new machinery is needed.

3. **Medium: s72 overwrites the primary EUR export artifacts.**

   Output paths contain chromosome, seed and release—but no cohort ([s72:91–94](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:91)). The outer loop consequently writes 3pop over EUR ([s72:185–186](/data/lvvyongji1/GWAS/code/s72_printed_matched.py:185)). The retained log confirms 720 samples, identifying the surviving 3pop export ([log:24–27](/data/lvvyongji1/GWAS/scratch/printed/s72_c22_0_0.log:24)).

   **Minimum fix:** include the cohort in export paths and preserve both cohorts’ files. The retained s68 EUR chr22 exports enabled my replay; the s72 artifact directory itself no longer preserves both cohorts.

4. **Lower severity: several affirmative statements still contradict the measurements or corrected interpretation.**

   “5/5 within-population TPR splits” is wrong ([claims:33](/data/lvvyongji1/GWAS/CLAIMS_FROM_RESULTS_ARB_20261017.md:33)): EAS is **4/5 on chr22** and **3/5 on chr1** ([raw chr22:2035](/data/lvvyongji1/GWAS/results/raw/s67_joint_passive_chr22.json:2035), [raw chr1:2034](/data/lvvyongji1/GWAS/results/raw/s67_joint_passive_chr1.json:2034)). Their mean intervals remain positive; this is a count correction.

   The off-causal reinterpretation is substantially fixed, but “sit at the nominal level” remains ([results:201–207](/data/lvvyongji1/GWAS/refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md:201)). That phrase still implies null calibration. “Individual QC is preserved” also exceeds the measured histogram and heterozygosity checks; partner heterozygosity changes remain in the JSON ([raw:38–41](/data/lvvyongji1/GWAS/results/raw/s71_protection_cost_chr22.json:38)).

   **Minimum fix:** correct the TPR counts, remove nominal-error comparisons, and state precisely which marginal quantities are preserved and which individual statistics change.

The real PLINK read-back and scatter/M normalization are correct. The reported s72 means and intervals also reconstruct correctly from its saved rows. The defect is the comparator construction and its interpretation—not substitution of float64 values for printed outputs.

### Memory Update

- **Closed:** per-split auditability, conditional-CI labeling, metric-specific s67 TPR oracle, and the principal J5 scope objection.
- **Open:** corrected printed-output superiority, particularly chr1 EUR; chr22 EUR’s repaired interval crosses zero.
- **Pattern:** calibration procedures and information grants drifted when the experiment moved into a new script.
- **Preserve:** the float64 result remains supported; positive-mean cross-population printed gaps remain inconclusive.