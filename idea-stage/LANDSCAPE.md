# Landscape Survey — Audit-Governed Multi-Site GWAS

**Direction:** Auditable reveal governance for MPC-based multi-site GWAS workflows — pre-registered reveal manifests, enforced reveal egress, cross-run release-relation auditing, empirical leakage measurement under named genomic attack models.

## 1. Sub-direction clusters

**A. Secure computation of GWAS statistics (the "how to compute" layer)**
- Ostrak, Randmets, Sokk, Laur, Kamm, "Implementing Privacy-Preserving Genotype Analysis with Consideration for Population Stratification," *Cryptography* 5(3):21, 2021. https://doi.org/10.3390/cryptography5030021 — MPC EIGENSTRAT/FastPCA/EMMAX/genomic control. Closest prior art to the proposal's execution layer.
- Aherrahrou, Tairi, Aherrahrou, "Genomic privacy preservation in GWAS: taxonomy, limitations, challenges, and vision," *Briefings in Bioinformatics*, 2024. https://doi.org/10.1093/bib/bbae356 — survey; catalogues MPC/HE/TEE/FL approaches and explicitly notes utility loss.
- PP-GWAS: "Privacy Preserving Multi-Site GWAS," *Nature Communications*, 2025. https://www.nature.com/articles/s41467-025-66771-z
- FedLMM/FedGLMM federated mixed-model GWAS, bioRxiv 2025. https://www.biorxiv.org/content/10.64898/2025.12.16.693409v1.full
- Ferreira et al., "A distributed computing solution for privacy-preserving GWAS (Gyosa)," bioRxiv 2024. https://doi.org/10.1101/2024.01.15.575678
- SECRET-GWAS: confidential computing for population-scale GWAS, bioRxiv 2024. https://doi.org/10.1101/2024.04.24.590989
- GA4GH TES federated meta-analysis: Chunduru et al., arXiv:2609.02227, 2026. [from proposal, unverified by this run]
- GDI D8.9 federated PoC, Zenodo, 2025. https://doi.org/10.5281/zenodo.16882880 [from proposal]

**B. Secure RELEASE of summary statistics (the "is the output safe" layer)**
- Pascoal, Decouchant, Völp, "Secure and distributed assessment of privacy-preserving GWAS releases (GenDPR)," *Middleware* 2022. https://doi.org/10.1145/3528535.3565253 — assesses whether a release is safe; the closest prior art to the proposal's release-risk evaluation.
- Pascoal PhD thesis (DyPS / I-GWAS / GenDPR), Univ. Luxembourg, 2022. https://orbilu.uni.lu/handle/10993/51996
- I-GWAS: "Privacy-preserving interdependent GWASes," *PoPETs* 2023(2). https://petsymposium.org/popets/2023/popets-2023-0026.pdf — explicitly: existing secure GWAS works "are insufficient to prevent privacy leaks once the GWAS result is made accessible"; adds interdependent membership checks across linked GWASes.
- "Enabling trade-offs in privacy and utility in genomic data sharing," *Genome Research*, 2023. https://genome.cshlp.org/content/early/2023/05/19/gr.277674.123.full.pdf — MI via Beacon responses and allele frequencies; query-suppression + noise as defense.
- "Privacy-Preserving and Efficient Verification of the Outcome in GWAS," PMC9536480 — quantifies MI risk from released statistics using Sankararaman LRT; compares risk vs metadata sharing.
- **PoPETs 2026 "Revisiting Assumptions for Membership Inference on Summary Statistics"** https://petsymposium.org/popets/2026/popets-2026-0134.php — shows attack accuracy reflects shared cohort/condition traits, not individual membership (≈10% separation, ~3x for TPR@1%FPR). Directly relevant: warns that measured MI signals can be confounded by group structure.

**C. Genomic membership/summary attacks (the measurement substrate)**
- Homer et al., *PLoS Genetics* 4(8):e1000167, 2008. Sankararaman et al., *Nature Genetics* 41(9):965–967, 2009. Shringarpure & Bustamante, *AJHG* 97(5):631–646, 2015. Raisaro et al., *JAMIA* 24(4):799–805, 2017. Bu et al., haplotype-based MI from summary stats, *Bioinformatics* 37(S1):i161, 2021. https://academic.oup.com/bioinformatics/article/37/Supplement_1/i161/6319711
- Wan et al., "Practical privacy protection of genomic data sharing through a differential privacy framework," *BMC Med Genomics* 2017 — per-individual privacy budget assigned in beacons.
- Humer / sequence-based MI survey work (recent; treat as [UNVERIFIED] for exact IDs).

**D. Language-level IFC / declassification / enforcement**
- Sabelfeld & Myers, "Language-Based Information-Flow Security," *IEEE JSAC* 21(1):5–19, 2003. Sabelfeld & Sands, "Dimensions and Principles of Declassification," *CSFW* 2005. Li & Zdancewic, "Downgrading policies and relaxed noninterference," *POPL* 2005.
- CISPA IFC group: "specification languages for information-flow policies" + "secure information-flow in web-based workflows" — runtime-monitoring privacy enforcement; policy-expressibility/decidability as the core challenge. https://cispa.de/en/research/research-areas/reliable-security-guarantees/information-flow-control
- Legacy-code IFC + declassification analysis (TU/e thesis). https://research.tue.nl/en/publications/information-flow-and-declassification-analysis-for-legacy-and-unt

**E. Query auditing / privacy odometers / budget**
- McSherry, "Privacy integrated queries (PINQ)," *SIGMOD* 2009. Roy et al., "Airavat," *NSDI* 2010. Rogers, Vadhan, Roth, Ullman, "Privacy Odometers and Filters: Pay-as-you-Go Composition," *NeurIPS* 2016.
- "GenoShare's Avatar Engine" — **simulatable auditing** for genomic requests: anticipate incoming queries, reply only to those that don't enable unauthorized inference; noted as limited to statistical (not individual-genotype) queries. (EPFL infoscience thesis chapter) https://infoscience.epfl.ch/bitstreams/2d0fd51e-8ab7-41d0-9652-91a8b0a3224c/download

**F. Workflow provenance / enforcement / attestation**
- W3C PROV-O (2013); Nextflow (Di Tommaso et al., *Nat Biotech* 2017); CWL; in-toto (USENIX Sec 2019); Certificate Transparency (RFC 6962). All **observe** ("prospective/retrospective provenance"); none **pre-authorize** an output before it happens.
- DataJoint + CWL: schema-driven provenance layer, "governance" but retrospective. https://finance.yahoo.com/sectors/technology/articles/datajoint-enables-seamless-migration-cwl-211000892.html
- SandboxEval (arXiv:2504.00018): taxonomy of untrusted-code sandbox escape/exfil behaviours (network, cloud upload, root, resource exhaustion) — useful as a violation-injection catalogue. https://arxiv.org/html/2504.00018v1
- LLM-assisted scientific workflow management (arXiv:2606.18425) — specification-driven workflow generation with explicit validation before execution; closest to "compile intent → validated artifact". https://arxiv.org/html/2606.18425v1

**G. Genomic privacy METRICS & their failures**
- Wagner, "Evaluating the Strength of Genomic Privacy Metrics," *ACM TOPS* 2016 — several widely used genomic privacy metrics are **non-monotonic** under adversary-strength increases ⇒ cannot anchor a security claim. https://doi.org/10.1145/3020003
- Ghasemian et al., "Safeguarding Privacy in Genome Research," PMC12150713, 2025 — bidirectional risk framework (summary stats + research outputs).

**H. Policy-as-code / DP pipeline enforcement**
- "Policy-as-code": governance rules executed inside pipelines, blocking non-compliant ops before production. Immutable/tamper-resistant audit trails are now standard data-governance practice (Atlan, GitLab pipeline-execution-policies). This is the *industrial* analogue of the proposal's enforced reveal egress, but at data-access granularity, not reveal-object-with-precision granularity.

## 2. What has been tried / what has NOT

**Tried:** secure computation of GWAS stats (A); assessing whether a *single* release is safe (B, DyPS/GenDPR/I-GWAS); genomic MI attacks and their DP beacon defenses (C); declassification theory and IFC type systems in-program (D); DP query accounting with noise (E); provenance logging (F); metric critique (G); industrial policy-as-code at access granularity (H).

**NOT tried (as a connected system):**
- Treating *every intermediate reveal of an MPC-GWAS workflow* (λ_GC, PCA spectrum, missingness, condition bits) as a first-class pre-registered object with precision/recipient/condition.
- **Enforcing** reveal egress at a narrow MPC-reveal choke point and measuring the enforcement boundary (blocked/detected/out-of-TCB) via adversarial injection.
- Detecting **cross-run release RELATIONS** (refinement, condition relaxation, receiver expansion, repetition, union-coverage) as a first-class audit object — as opposed to per-run safety checks or per-individual DP budgets.
- Binding the above to **genomic attack-model leakage estimates** as an empirical frontier, honestly separated from the syntactic (provable) layer.

## 3. Recurring Future-Work limitations
- "Privacy-utility trade-off" stated as unsolved in every FL/DP GWAS survey; noise destroys GWAS power.
- DyPS/GenDPR family: safe **single** release; no treatment of workflow-level intermediate reveals or cumulative multi-run release.
- IFC: enforcement is in-program; no cross-organizational workflow object with domain semantics.
- Provenance: explicitly observational; cannot refuse an unauthorized output.
- Beacon/DP budgets: per-individual ε accounting; assumes an established DP sensitivity, which is unavailable/infeasible for high-dimensional sparse genotypes with GWAS-level utility.

## 4. Open problems named by multiple sources
1. Post-computation leakage is under-addressed: "insufficient to prevent leaks once the result is accessible" (I-GWAS, GenDPR, PP-GWAS surveys).
2. No principled utility-preserving enforcement for GWAS summary release (all surveys).
3. Metrics don't compose / aren't monotone (Wagner; PoPETs 2026 MI-assumption paper).
4. Repeated/interdependent releases compound (I-GWAS "interdependent"; beacon repeated-query attacks; DP composition) — but no workflow-level cross-run *relation* detector.

## 5. Structural gaps (candidate-generation seeds)

- **G1 (method-transfer):** Declassification-policy enforcement (D) is untried as a *workflow-level reveal-object* governance layer in genomic multi-site analysis.
- **G2 (method-transfer):** Transparency-log integrity (CT/in-toto, F) is borrowed for tamper-evidence but never to prove an output was *pre-authorized before it occurred*; the "authorization-before-event" property is uninstantiated for reveal objects.
- **G3 (untested-assumption):** "Intermediate diagnostic reveals (λ_GC, PCA spectrum, missingness) have negligible marginal leakage beyond final summary stats" is assumed, never measured under a named attack model.
- **G4 (untested-assumption):** "MPC's reveal/open is the uniquely narrow enforcement point" — assumed, not compared against logging-only provenance and federated-instrumented baselines under identical bypass injection.
- **G5 (contradiction):** Wagner's non-monotonicity (G) vs the implicit "finer precision ⇒ more leakage" design assumption. Do precision/condition changes on genomic reveals behave monotonically? The proposal bets both ways and hedges.
- **G6 (diagnostic):** Do cross-run release *relations* (refinement / condition relaxation / receiver expansion / union-coverage) track attack-risk growth at all — or is the ledger just provenance with extra vocabulary?
- **G7 (untested-assumption):** Pre-commitment (moving judgment from run-time to compile-time) has an acceptable GWAS statistical cost — sign unspecified.
- **G8 (diagnostic/contradiction):** PoPETs 2026 says measured MI "gains" often reflect shared cohort traits, not membership. Does the proposal's leakage frontier survive ancestor-matched, panel-disjoint, condition-controlled evaluation?

## 6. Closest existing claim to the proposal's combination

**None found that combines all four** (pre-registered reveal manifest + enforced egress + cross-run release-relation auditing + genomic attack evaluation).
- Closest single-system: **GenDPR** (release safety assessment) and **I-GWAS** (interdependent releases) — but both assess/limit *which SNPs are released*, not the *workflow's full reveal set with enforced pre-registration and cross-run relation audit*.
- Closest mechanism: **declassification allowlist + taint/IFC** (D) and **policy-as-code + immutable audit** (H) — but neither carries genomic reveal-object semantics nor attack-model binding.
- Closest measurement: genomic MI attacks (C) + Wagner metric critique (G).

⇒ The proposal's novelty is *integration + honest boundary measurement*, not a new cryptographic or privacy-definitional primitive; this is the main novelty threat to pre-empt (see proposal §9 Risk 2).

## 7. Verification status
Papers marked with a URL were returned by live search (Tavily) in this run. Entries tagged [UNVERIFIED] or marked "[from proposal]" were not independently re-verified in this run and MUST be checked with `verify_papers.py` before citation in the paper.
