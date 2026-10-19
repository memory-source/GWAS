# Idea Discovery Report — Supplemental Search Round (Novelty Re-confirmation)

**Direction**: Re-run the adjacent-work search around the multi-site-GWAS reveal-governance
paper (PC-loadings as a measured intermediate leak; governed release manifest / egress
monitor / hash-chained ledger; cross-run composition; conditioned-AF vs PC ordering).
**Date**: 2026-10-13
**Scope**: 2024–2026 adjacent work that could displace or demand positive handling of the
paper's novelty statements.
**Pipeline (lite)**: research-lit → adjacency triage → novelty re-check. Idea generation /
pilots do not apply: the paper and its codebase already exist; this round only re-tests the
novelty claim, not a new idea.

## Executive Summary
The paper's novelty statement survives the supplemental search. Five families were
re-touched: (1) GWAS membership inference from summary statistics, (2) released-PCA /
loadings leakage, (3) multi-site / federated GWAS secure computation and its released
objects, (4) release provenance / egress control, (5) benchmark-adjacent work. One real gap
was found and fixed in the paper: the multi-site secure-GWAS literature (SF-GWAS,
*Nat Genet* 2025) shares study-fitted PC loadings across sites, which is exactly the object
this paper measures as the stronger intermediate leak; it was not cited. Two further
adjacent leakage results (released-PCA reconstruction, GWAS-summary genotype recovery) were
also uncited. All three are now cited and positively scoped in `paper/main.tex`
§Related Work; no other entry changes a headline claim.

## Literature Landscape (adjacent families)

### 1. GWAS membership inference from summary statistics
- **Berrang et al., PoPETs 2026** (already cited `popets2026`): reported MI gains largely
  reflect shared cohort traits, not membership. The paper's disjoint-negative control
  directly addresses this; unchanged.
- **Lan, Pawitan, Shen (medRxiv 2025.09.09)**, "Genomic privacy risks in GWAS summary
  statistics" — linear-programming genotype recovery from summary statistics plus
  high-dimensional phenotype data; full recovery at phenotype/sample ratio > 0.85,
  re-identification above 0.16; non-European populations more vulnerable. **Not a
  membership-inference competitor** (it recovers genotypes, a strictly stronger objective,
  under an additional phenotype-data assumption) and **not a governed-release competitor**
  (it treats the summary-statistics release, not the workflow's intermediate reveal set).
  NEW CITATION `lan2025`; positive scope in §Related Work.

### 2. Released-PCA / loadings leakage
- **Zari et al., PSD 2022** (already cited `zari2022`): MIA against released principal
  components; stronger for small training sets and more components. This is the source
  primitive the paper reuses; unchanged.
- **Kwatra & Torra, SocialSec 2023**: data-reconstruction attack against PCA, extending
  Zari's MIA to a strictly stronger objective (full data reconstruction), incl. from
  synthetic-data-derived components, with DPPCA as defense. **No competitor**: different
  objective, no GWAS setting, no governance. NEW CITATION `kwatra2023`; strengthens the
  paper's "not a new attack primitive" scoping.
- **Genomic PCA / PCA-QC tooling** (genomicPCA, GCTA/smartpca/FlashPCA docs): confirms the
  loadings are routinely computed and shared; background, not a competitor.

### 3. Multi-site / federated GWAS secure computation
- **SF-GWAS, Froelicher et al., *Nature Genetics* 2025**: secure federated PCA- and
  LMM-based GWAS at biobank scale via MPC+HE; each site runs the study-fitted PCA and
  **shares the loadings (and small intermediates) across sites** so a site can project its
  own samples. This is the same "federated PCA shares the loadings" release the paper
  measures as its stronger intermediate channel; it was previously asserted in the paper
  without a citation. NEW CITATION `sfgwas2025`; the paper now cites it as the concrete
  release object it measures. **No competitor**: SF-GWAS secures the computation and does
  not measure the leakage of the released diagnostic set, which is the paper's object.
- **PP-GWAS (Nat Commun 2025)**, **FedLMM (2025)**, **Ostrak 2021**, **Aherrahrou 2024
  survey** (all already cited): secure-computation layer; unchanged.

### 4. Release provenance / egress control
- **TRE / egress-airlock practice** (Lifebit, DNAnexus TRE AUP 2026, AWS DRTC, Five Safes /
  SATRE 2.0): production genomic platforms enforce controlled egress and keep governance
  lineage. These are **engineering controls, not measured governance research**: none binds
  a pre-registered reveal manifest to an attack-model evaluation or audits cross-run
  release-relation composition. The paper's contribution sits between provenance logging
  (`prov-o`, `intoto`, `rfc6962`, all cited) and these production controls. No new citation
  required; the paper's Related Work already names the logging side.

### 5. Benchmark-adjacent
- LLM-agent orchestration benchmarks (BioAgent Bench 2026, BioML-bench, LAB-Bench, the
  2026 LLM-agent-evaluation survey) evaluate *task completion / pipeline building*. None
  measures the leakage of a pipeline's intermediate diagnostic **release set** under a
  membership adversary; they are orthogonal evaluation substrate. No citation needed.

## Novelty Verification (re-confirmation)
**Claim under test**: "a controlled measurement of the *intermediate* releases (PC loadings,
λ_GC, per-site AF, missingness, haplotype) that the secure-computation and release-safety
literatures leave unmeasured, plus a governance layer over the release set."

- **Closest system**: SF-GWAS (shares the loadings; does not measure their leak) and
  GenDPR/I-GWAS (assess *which SNPs* a single release may expose; not the workflow's full
  intermediate reveal set). **Differentiation holds**: the paper's object is the released
  diagnostic objects and their *ordering*; no located work measures the PC-vs-AF ordering
  on the released set, nor the cross-run composition of manifest-permitted releases.
- **Closest attack primitive**: Zari 2022 (PCA MIA) and Kwatra 2023 (PCA reconstruction).
  The paper does not claim the primitive; it reuses it in a controlled comparison. Scoping
  in §Related Work is explicit and now cites both.
- **Closest metric critique**: Wagner 2016 (already cited) and Berrang 2026 (already
  cited). Unchanged.
- **Verdict**: **PROCEED** — the novelty statement is a *comparison + governance* claim,
  not a primitive claim, and the search found no work that occupies that comparison. The
  one actionable finding (SF-GWAS uncited) is fixed.

## Changes applied to the paper
- `paper/main.tex` §Related Work: added SF-GWAS as the concrete federated-loadings-sharing
  system the paper measures, plus the released-PCA-reconstruction and GWAS-summary
  genotype-recovery neighbors, each with an explicit "not our contribution" scoping clause.
- Added bibliography entries `sfgwas2025`, `kwatra2023`, `lan2025`.
- Recompiled: 21 pp, 0 errors, 0 undefined references.
- No headline number, claim, or scoping sentence outside §Related Work changed.

## Eliminated / non-adjacent
- LLM-agent orchestration benchmarks (task-completion focus) — orthogonal.
- Production TRE/egress-platform marketing docs — engineering controls, not measured
  governance research; the logging side is already cited.
- Genomic PCA-QC tooling docs — background.

## Next Steps
- [ ] No further experiments from this round; it re-confirms novelty and adds three
      citations. Existing claim/limitation closure (round 16) stands.
