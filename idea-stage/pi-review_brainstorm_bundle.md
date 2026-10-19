# Idea-generation bundle — audit-governed multi-site GWAS (pi-review brainstorm)

You are a senior ML/security researcher brainstorming research ideas. This is a
**generator** task: produce a broad, bold pool of candidate ideas. Do NOT rank or
filter — strict filtering happens in a separate triage stage.

## Research direction
Auditable reveal governance for MPC-based multi-site GWAS workflows:
pre-registered reveal manifests, enforced reveal egress at the MPC reveal
primitive, cross-run release-relation auditing, and empirical leakage measurement
under named genomic attack models.

## Current landscape (live-search verified unless tagged)
- **Secure-compute layer (crowded):** Ostrak et al. 2021 (MPC EIGENSTRAT/EMMAX); Aherrahrou et al. 2024 survey (BiB); PP-GWAS (Nat Commun 2025); FedLMM/FedGLMM (bioRxiv 2025); Gyosa (bioRxiv 2024); SECRET-GWAS (bioRxiv 2024); GA4GH TES federated meta (arXiv 2026); GDI PoC (2025). These all stop at "compute the statistic securely".
- **Release-safety layer:** Pascoal/Decouchant/Völp GenDPR (Middleware 2022) assesses whether a *single* release is safe; I-GWAS (PoPETs 2023) does interdependent membership checks; Genome Research 2023 does Beacon/AF query-suppression + noise; PoPETs 2026 "Revisiting Assumptions for MI on Summary Statistics" shows measured MI gains often reflect shared cohort traits not individual membership.
- **Attack substrate:** Homer 2008; Sankararaman 2009; Shringarpure & Bustamante 2015; Raisaro 2017; Bu 2021 (haplotype MI); Wan 2017 (DP beacon budget).
- **IFC/declassification:** Sabelfeld & Myers 2003; Sabelfeld & Sands 2005 (declassification dimensions); Li & Zdancewic 2005; CISPA runtime-monitoring IFC; GenoShare "simulatable auditing" (limited to statistical, not genotype, queries).
- **Query auditing:** PINQ 2009; Airavat 2010; privacy odometers (Rogers et al. 2016). All assume a DP sensitivity/budget; no such budget exists for sparse high-dim genotypes with GWAS utility.
- **Provenance/enforcement:** W3C PROV; Nextflow/CWL; in-toto; Certificate Transparency (RFC 6962); SandboxEval (arXiv:2504.00018) taxonomy of sandbox escapes; LLM workflow-management (arXiv:2606.18425) compiles intent→validated artifact. All provenance is *observational* (log after the fact); none *pre-authorize* an output.
- **Metric critique:** Wagner 2016 (ACM TOPS) — common genomic privacy metrics are non-monotone under adversary strength.

## Key structural gaps
- G1 declassification enforcement untried as workflow-level *reveal-object* governance in genomics.
- G2 transparency-log integrity borrowed, but "output was pre-authorized before it occurred" is uninstantiated for reveal objects.
- G3 "intermediate diagnostic reveals (λ_GC, PCA spectrum, missingness) have negligible marginal leakage" — assumed, never measured under a named attack.
- G4 "MPC reveal is uniquely narrow enforcement point" — assumed, never compared to logging-only provenance + instrumented federated baselines under identical bypass injection.
- G5 Wagner non-monotonicity vs the "finer precision ⇒ more leakage" design assumption.
- G6 do cross-run release *relations* (refinement/condition-relax/receiver-expand/union-coverage) track attack-risk growth at all?
- G7 pre-commitment's GWAS statistical cost — sign unknown.
- G8 does the leakage frontier survive ancestor-matched, panel-disjoint, condition-controlled evaluation (per PoPETs 2026)?

## Known novelty threat
All four ingredients trace to known work; a systems reviewer will call it "policy
plumbing" unless there is a falsifiable measurement, a formal statement, or an
instantiated property that is genuinely new. The direction's bet is on **honest
boundary measurement + integration**, not a new primitive.

## Your task
Generate **8-12 concrete research ideas**. For each:
1. One-sentence summary
2. Core hypothesis (what you expect and why)
3. Minimum viable experiment (cheapest way to test — this is a CPU/MPC project, NOT GPU)
4. Expected contribution type: empirical finding / new method / theoretical result / diagnostic
5. Risk level: LOW / MEDIUM / HIGH
6. Estimated effort: days / weeks / months

Prioritize ideas that are:
- Testable with moderate compute (≥16-core CPU, no large GPU)
- Likely to produce a clear positive OR negative result (both publishable)
- Simple at the core: one mechanism, few moving parts. If novelty only appears once a second module is added, that is packaging, not novelty.
- Aware of the landscape above — awareness, not avoidance. Differentiation is the novelty check's job later.

"Apply X to Y" is legitimate when the application reveals something non-obvious.
Do not steer around crowded areas — proximity to strong work means the problem
matters. Be genuinely creative: surprising connections, inverted assumptions,
questions nobody thought to ask. Generate first, filter later.

Return the ideas as a numbered markdown list. Do not write files.
