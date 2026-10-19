# Research Idea Report — Audit-Governed Multi-Site GWAS

**Direction**: Auditable reveal governance for MPC-based multi-site GWAS workflows — pre-registered reveal manifests, enforced reveal egress, cross-run release-relation auditing, empirical leakage measurement under named genomic attack models.
**Base document**: `proposal_experiment_plan.md` (full PoPETs-style systems + empirical study).
**Generated**: 2026-10-02
**Ideas evaluated**: 17 generated (12 cross-model brainstorm on glm-flash + 5 executor-lens extras) → 17 passed the feasibility gate (all CPU-only) → 3 selected for the pilot slate → 3 recommended.
**Validation mode**: paper validation only (no GPU); the "pilots" below are CPU feasibility decisions, not training runs.

## Landscape Summary

The multi-site-GWAS field has a mature, crowded **secure-computation** layer: MPC/HE/TEE pipelines that compute GWAS statistics without pooling genotypes (Ostrak et al. 2021; Aherrahrou et al. 2024 survey; PP-GWAS, *Nat Commun* 2025; FedLMM, 2025; Gyosa, 2024; SECRET-GWAS, 2024; GA4GH-TES federated meta, 2026). All of it stops at "compute the statistic securely."

A second layer asks whether a **single release** is safe: GenDPR (Pascoal et al., *Middleware* 2022) assesses GWAS release privacy; I-GWAS (*PoPETs* 2023) adds interdependent membership checks across linked GWASes and explicitly notes that secure-GWAS work "is insufficient to prevent privacy leaks once the result is made accessible"; the *Genome Research* 2023 line adds query-suppression + noise. The **attack substrate** is well established (Homer 2008; Sankararaman 2009; Shringarpure & Bustamante 2015; Raisaro 2017; Bu 2021), and Walters/2026 PoPETs work warns that measured MI "gains" are frequently confounded with shared cohort/ancestry traits rather than individual membership.

Two adjacent mature traditions supply the mechanism the proposal wants to borrow: **language-level IFC / declassification** (Sabelfeld & Myers 2003; Sabelfeld & Sands 2005; Li & Zdancewic 2005; CISPA runtime-monitoring IFC; GenoShare's "simulatable auditing") and **query auditing / provenance / attestation** (PINQ 2009; Airavat 2010; privacy odometers 2016; W3C PROV; Nextflow/CWL; in-toto; Certificate Transparency RFC 6962). Crucially, all provenance/attestation systems are **observational** — they log after the fact and cannot refuse an unauthorized output — and all query-auditing systems assume a DP sensitivity that arguably does not exist for sparse, high-dimensional genotypes at GWAS utility.

The gap the proposal occupies: nobody has connected *pre-registered reveal manifests + enforced reveal egress at a narrow choke point + cross-run release-relation auditing + genomic attack-model leakage measurement* into one auditable system, and—more importantly—nobody has **measured the boundary** of any of those four claims. Wagner (2016) shows the field's own privacy metrics are non-monotone under adversary strength, which forbids the "budget" framing and forces the honest syntactic-vs-empirical split the proposal already adopts. The novelty threat is real and named: a systems reviewer will call the mechanism "a domain-specialized declassification allowlist plus a transparency log," so the paper must live on honest boundary measurement and integration, not on a new primitive.

## Recommended Ideas (ranked)

### Idea 1: Matched-control leakage frontier for GWAS reveal objects (T5 — credibility gate)
- **Method (what we actually do)**: (1) Build a multi-site panel from 1000 Genomes (and msprime as a second, controllable source), split into pseudo-sites with known ancestry composition, run a standard PLINK2 QC → PCA → association pipeline. (2) Implement two membership-inference attack families (an aggregate allele-frequency Homer/Sankararaman likelihood-ratio attack; a Beacon/AAF-style likelihood-ratio or small-ML attack). (3) Run each attack on each release object (F1–F5 × M0–M3), but with hard controls: members and non-members are ancestry-matched, the attacker's reference panel is disjoint from the target panel, LD-pruned sites, and an explicit "ancestry-only" control attack (label/PCA-projection) as an upper bound on the confound. (4) Decompose measured AUROC / TPR@FPR into an ancestry-signal component (stratified permutation) vs an individual-membership component.
- **Hypothesis**: A large fraction of reported summary-stat leakage is ancestry/cohort confound; the *residual* individual-membership leakage concentrates in the low-MAF tail and in effect-size (not p-value) vectors.
- **Minimum experiment**: 1000G-derived 3–5-site split, ~2k total individuals, chr22 or ~10k LD-pruned sites, both attack families, ancestry-matched + panel-disjoint controls, stratified-permutation decomposition.
- **Expected outcome**: Either (a) residual leakage survives controls and is localized to specific objects/precision → the frontier is credible and the governance story holds; or (b) residual leakage vanishes under control → an honest, citable calibration result that recalibrates the whole attack literature and still supports the paper's "empirical layer is separate from syntactic layer" thesis.
- **Novelty**: 8/10 — closest work: PoPETs 2026 "Revisiting Assumptions for MI on Summary Statistics" warns of the confound; **nobody has run the controlled decomposition for GWAS *release objects* across a release-policy grid.** Differentiable.
- **Feasibility**: CPU-only, ≥16 cores, ~64 GB RAM; 1000G chr22 subset + PLINK2 + scikit-allel + statsmodels; weeks. No GPU.
- **Risk**: LOW
- **Contribution type**: empirical finding
- **Pilot result**: PAPER-VALIDATED: mandatory precondition for every leakage claim; if this fails, all privacy conclusions in the proposal are blocked (proposal P0-3 fail branch). Decision: **must run first.**
- **Reviewer's likely objection**: "This is calibration infrastructure, not a contribution." — Answer: it is the gate that makes every downstream number citable, and a *clean negative* here is itself a field-level result.
- **Why we should do this**: Nothing else in the project is believable until the leakage meter is trustworthy; it is low-risk and shares its entire build with Ideas 2 and 3.

### Idea 2: Marginal leakage of intermediate QC diagnostics (T2 — the paper's spine)
- **Method (what we actually do)**: (1) Reuse Idea 1's panels and attacks. (2) Fix a final-release policy (conservative F1/F2 and permissive F3/F4) and, one at a time, add each intermediate diagnostic reveal that real multi-site GWAS pipelines emit: λ_GC alone; PCA eigenvalue spectrum + PC–phenotype correlations; per-SNP missingness / QC tables; and the secret-dependent branch *condition bit* (which PCA count was chosen). (3) Measure the attack-success *increment* over the previous tier, under the ancestry-matched / panel-disjoint controls from Idea 1, reporting TPR@FPR≤0.1% and ≤1% per ancestor.
- **Hypothesis**: PCA spectrum/loadings alone support above-chance site membership inference (eigen-decomposition encodes kinship/ancestry structure), while λ_GC alone leaks little; the "only publish final p-values is low-leakage" intuition is incomplete because the middle reveals carry independent signal.
- **Minimum experiment**: 2 final-release tiers (conservative + permissive) × 4 intermediate tiers (M0–M3) + condition-bit on/off, using Idea-1 attacks.
- **Expected outcome**: Positive (some diagnostic has a direction-consistent, measurable marginal effect) → "intermediate reveals are a real budget" and the governance mechanism is motivated. Negative (final release dominates; marginal effect negligible) → reframe the paper as **"an empirical release-strategy frontier for multi-site GWAS"**, keeping intermediate reveals as a documented negative result (proposal P0-4 fail branch). Both are publishable.
- **Novelty**: 8/10 — GenDPR/I-GWAS assess *final* releases only; marginal leakage of *intermediate* workflow diagnostics is unmeasured. Differentiable.
- **Feasibility**: CPU-only, shares Idea-1 build; weeks.
- **Risk**: LOW–MEDIUM (the effect may be ancestry-confounded — but that is exactly what Idea 1's controls rule out).
- **Contribution type**: empirical finding
- **Pilot result**: PAPER-VALIDATED: this is proposal contribution **C3**, the paper's center of gravity; it decides the narrative (positive vs reframe). Decision: **run immediately after Idea 1.**
- **Reviewer's likely objection**: "Any PCA-spectrum effect is just ancestry re-detection." — Named risk, answered only by surviving Idea-1's matched controls; not patched with an extra module.
- **Why we should do this**: It is the single highest-information experiment — one build, two publishable outcomes, and it sets the whole paper's direction.

### Idea 3: The sanitized-diagnostics paradox (T11 — the high-upside bet)
- **Method (what we actually do)**: (1) Define a scripted "analyst" that performs GWAS QC under a fixed reveal budget and, when it needs a diagnostic it does not have, issues *probe* reveals (thresholded re-queries, finer-precision re-releases) to reconstruct it. (2) Run two independent scripted policy families × ≥2 panel seeds, each with and without *sanitized* diagnostics pre-provided (binned PCA spectrum, rounded λ_GC, coarsened missingness). (3) Count total probe reveals and cumulative measured attack success (using Idea 1's meter) in each arm, at equal analysis quality (same final hits / power).
- **Hypothesis**: Pre-providing sanitized diagnostics *reduces* cumulative risk by removing the analyst's incentive to probe — the total-risk curves cross, i.e., governance should sometimes *add* a well-chosen reveal rather than suppress everything.
- **Minimum experiment**: pre-registered comparison, 2 policy families × 2 seeds × 2 arms; primary endpoint = probe-reveal count; secondary = cumulative attack success at equal analysis quality.
- **Expected outcome**: Positive if both policy families show ≥20% fewer probe reveals at equal analysis quality, consistent across seeds → a counterintuitive, memorable result with zero prior work. Negative if the effect is policy-family-dependent → dropped as a policy artifact (recorded, not rescued by tuning policies).
- **Novelty**: 10/10 — no prior work at all; an inverted framing of the standard "suppress to protect" instinct.
- **Feasibility**: CPU-only; reuses Idea-1 leakage meter; scripted policies in Python; weeks.
- **Risk**: HIGH (simulation-validity threat: the analyst model is scripted, so a reviewer will ask whether the probing behavior is realistic).
- **Contribution type**: empirical finding (counterintuitive / diagnostic)
- **Pilot result**: PAPER-VALIDATED with a **named conviction criterion**: dropped unless ≥20% probe-reveal reduction at equal quality appears in *both* independent policy families and ≥2 seeds.
- **Reviewer's likely objection**: "The probing behavior is an artifact of your scripted policy." — This is the idea's named risk; the mitigation is two independent policy families and a pre-registered drop rule, *not* adding realism modules.
- **Why we should do this**: The pool's highest raw upside and its designated non-safe pick; if it holds it is the paper's most memorable finding, and it is cheap enough to try after the spine is secured.

## Selected-for-pilot slate (per cross-model triage)
Rank order: **1. T5 (matched-control frontier)** → **2. T2 (marginal intermediate leakage)** → **3. T11 (sanitized-diagnostics paradox)**.
T5+T2 share one build and match the proposal's own gating (P0-3 → P0-4); T11 is the discretionary high-upside slot. The systems bake-off items (T15 enforce-vs-detect-vs-noisy, T7 bypass-injection) are **not** pilot slots — they are executed inside the already-committed P0-5 build once the measurement spine is settled.

## Eliminated / deferred ideas (for reference)
| Idea | Status | Reason |
|------|--------|--------|
| T10 composition-additivity theorem | deferred | worst effort/risk ratio for this budget; fixed-panel case is corollary-shaped vs DP folklore; adaptive violation map is simulation |
| T8 slicing attack on gate granularity | deferred | likely **out of the A2 threat model** (constrained authors may never be authorized enough scalar reveals); risk of vacuous null — revisit after T15 fixes the threat model |
| T13 TEE-vs-MPC parity | deferred | low information yield; simulated SGX cannot test the side-channel trust story; proposal already keeps TEE as a qualitative table |
| T12 receiver-dimension collapse | cheap add-on | shallow but cheap; keep as a supplement audit if time permits |
| T1 DP-inadmissibility map | supporting material | tutorial-shaped; keep as a table justifying non-DP governance, not a headline |
| T6 manifest grammar coverage | supporting | design-defining but corpus-bias risk; feeds the DSL spec |
| T3/T17 precision & condition-bit | fold into T2 | mandatory anyway (P0-8 sub-items); run inside Idea 2's grid |
| T4 relations-track-risk | gated on T2 | run only if T2 shows signal; otherwise C4 demoted per P0-7 |
| T7/T15 enforce-vs-detect-vs-log | in P0-5 | highest-persuasion systems panel; no separate pilot slot needed |
| T9 cost of pre-commitment | P1 | important but P1-level; null outcome weakly informative |
| T14 static decidability | supporting | formal core for the DSL checker; moderate surprise |

## Suggested execution order
1. **Idea 1 (T5)** — build panels + attack implementations + matched-control decomposition. Everything depends on this meter.
2. **Idea 2 (T2)** — add the intermediate-reveal grid on top of Idea-1's build; decides positive vs reframe.
3. **Idea 3 (T11)** — scripted-analyst experiment with the pre-registered drop rule.
4. In parallel (already committed): P0-5 enforce/detect/log bake-off (T7/T15), then P0-7 relation-risk regression (T4) only if T2 has signal.

## Next Steps
- [ ] Freeze the measurement spine (panels, attacks, controls) and lock it in `idea-stage/RESEARCH_CONTRACT.md`.
- [ ] Invoke `/experiment-bridge` to implement Idea 1 → Idea 2 → Idea 3 (see the plan below).
- [ ] After first results, run `/result-to-claim` before writing any claim.
- [ ] Write paper draft, then `/auto-review-loop`.

## File management (this stage)
- `idea-stage/IDEA_REPORT.md` — this file (canonical deliverable of `/idea-creator`).
- `idea-stage/IDEA_CANDIDATES.md` — merged, annotated candidate set (T1–T17).
- `idea-stage/LANDSCAPE.md` — Phase-1 landscape + structural gaps.
- `idea-stage/pi-review_brainstorm_bundle.md`, `pi-review_triage_bundle.md` — cross-model input bundles.
- `.aris/traces/idea-creator/2026-10-02_run01/` — full prompt/response traces (3 pi-review calls).
- Downstream stages write to `refine-logs/`, `review-stage/`, `paper/`, `code/`, `results/` per the Output Versioning Protocol.
