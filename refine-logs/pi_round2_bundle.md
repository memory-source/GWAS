# Cross-model re-review bundle — Audit-Governed Multi-Site GWAS (round 2)

You reviewed the round-1 bundle and returned C1 PARTIAL / C2 PARTIAL / C3 NOT
SUPPORTED. The executor fixed the flagged issues and re-ran. Judge the REVISED
evidence. Full detail in refine-logs/EXPERIMENT_RESULTS.md and results/raw/.

## What changed since round 1
1. **Disjoint negatives.** All negative sets are now an INDEPENDENT non-member
   panel disjoint from BOTH the study and the attacker's reference panel. This
   removed the self-referential bias you flagged.
2. **Bootstrap CIs** added (95% CI on TPR@FPR≤1%).
3. **S3 redesigned** with a non-degenerate coarse-first search, a measured
   quality function, and a documented scripted-analyst model.
4. F1 (degenerate) tier dropped; a real site-size × precision boundary curve
   (S5) added.
5. msprime cross-check rerun.

## Revised results
**S1b frontier (disjoint negatives, 2 seeds):** AUROC — F2 top-200 0.578;
F4 AF-1dp 0.517; F3 AF-2dp 0.712; F5 full 0.720; F5 per-site 0.859.
TPR@1%: 0.024 / 0.014 / 0.061 / 0.065 / 0.177 with 95% CIs (e.g. F5
[0.041,0.113]). Ancestry-only control 0.495 at all tiers. **Fixing the negative
set roughly halved the earlier AUROCs (F5 0.92→0.72).**

**S5 boundary (TPR@1%):** n=50: 0.09/0.61/0.60 (1dp/2dp/full); n=100:
0.03/0.38/0.36; n=200: 0.02/0.11/0.10; n=400: 0.02/0.03/0.03; n=800:
0.01/0.03/0.03.

**S2b intermediate reveals (disjoint negatives, 3 seeds, 5-fold CV):** final-only
0.753; +λ_GC 0.753; +condition bit 0.753; +missingness 0.754; **+PCA spectrum
0.861 (Δ +0.108)**.

**S3 paradox (scripted analyst):** probe reveals 3→0 (family A) and 4→1
(family B) with sanitized diagnostics, at equal quality.

**S4 msprime:** full-AF AUROC 0.716, AF-1dp 0.627.

## Your task
1. Re-verdict C1/C2/C3 as SUPPORTED / PARTIAL / NOT SUPPORTED.
2. Did the disjoint-negative fix change your confidence, and in which direction?
3. Is Δ+0.108 (PCA spectrum) now credible, or is a residual artifact present?
4. Rank the claims by value to a reviewer.
5. Single most important remaining experiment before submission; and the
   strongest surviving rejection argument (do not propose adding modules).
