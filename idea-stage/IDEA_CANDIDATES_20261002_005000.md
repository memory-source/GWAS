# Idea candidates — merged, deduped, annotated (Phase 3)

Sources merged: pi-review brainstorm (12 candidates, model glm-flash, thread
673924b852ba46c5b27c0326eafb980b) + executor lens fan-out (Tier-3 inline, 6
candidates). Mechanical dedup only — no quality filtering. Every candidate is
automatically project-feasible for the feasibility gate: all are CPU-only
(≥16-core CPU; msprime/PLINK2/SPy available), none needs a GPU or a
provably-unavailable dataset.

Legend: T# = merged id. `prior_work` = nearest existing work (annotation for the
jury, not a filter). `so_what` = why the result matters either way.

---

## T1 — DP-inadmissibility map of GWAS reveal artifacts
- Hypothesis: standard reveals split into bounded-sensitivity (λ_GC, missingness) vs unbounded (PCA loadings, GRM/kinship, raw p-vectors), so no DP budget admits the latter at GWAS utility ⇒ non-DP egress governance is the only remaining control point.
- MVE: closed-form + numerical substitution-sensitivity per artifact type on synthetic panels (n≈1k, m≈100k). Effort: days–weeks. Risk: LOW.
- prior_work: PINQ/Airavat assume DP sensitivity; Wagner 2016 (metric critique). No one has published the artifact-by-artifact sensitivity dichotomy for GWAS reveals.
- so_what: even a purely negative/confirmatory map redirects the field — it tells governance designers *which* reveals a DP budget could ever cover.
- Effort_note: pure analysis, no MPC.

## T2 — Marginal leakage of intermediate QC diagnostics  ★ (proposal C3 / P0-4)
- Hypothesis: PCA eigen-spectrum/loadings alone support above-chance site membership inference; λ_GC alone leaks little. Replaces an assumed "negligible" with a measured per-artifact boundary.
- MVE: 1000G split into pseudo-sites; standard QC; feed each diagnostic subset to Sankararaman/Homer attacks vs full-summary ceiling; AUROC per artifact. Effort: weeks. Risk: LOW–MEDIUM.
- prior_work: GenDPR/I-GWAS assess release safety of *final* releases only; nobody measures marginal leakage of *intermediate* diagnostics. PoPETs 2026 warns MI gains can be cohort-confounded.
- so_what: decides the paper's spine. Positive ⇒ "middle reveals are a real budget"; negative ⇒ reframe as a publication-strategy frontier. Either outcome is publishable.

## T3 — Precision non-monotonicity of reveal objects (proposal P0-8 / G5)
- Hypothesis: coarser precision sometimes *increases* attack success (censoring plateaus become fingerprintable), contradicting "round more ⇒ leak less".
- MVE: precision grid (5–8 levels × 3 artifact types) × 2 attacks on fixed panels; AUROC vs precision. Effort: days–weeks. Risk: MEDIUM.
- prior_work: Wagner 2016 non-monotonicity of *metrics*; Sanders/others in DP. Untested for GWAS reveal-object precision as a *control knob*.
- so_what: if monotone, precision becomes a defensible knob; if not, the proposal's precision lattice must be demoted to audit metadata — both are clean findings.

## T4 — Do cross-run release relations track attack-risk growth? (proposal P0-7/C4 / G6)
- Hypothesis: union-coverage and refinement edges predict measured AUROC growth; receiver-expand alone does not (receivers don't change the statistic).
- MVE: ~50–100 ordered release pairs; label relation; regress AUROC delta on relation indicators. Effort: weeks. Risk: MEDIUM.
- prior_work: I-GWAS "interdependent" membership checks; DP composition; no workflow-level *relation*-to-risk mapping.
- so_what: tells whether the ledger is provenance-with-vocabulary or genuinely risk-bearing. Decides whether C4 stays a main contribution.

## T5 — Leakage frontier under matched controls (proposal P0-3 / G8)
- Hypothesis: much reported summary-stat leakage is ancestry/cohort confound; residual concentrates in the low-MAF tail and effect-size (not p-value) vectors.
- MVE: multi-site 1000G simulation, ancestry strata + matched controls; Homer / Sankararaman / LRT; decompose AUROC into ancestry vs individual components by stratified permutation. Effort: weeks. Risk: LOW.
- prior_work: PoPETs 2026 (MI reflects shared cohort traits); Venkatesaramani beacon defenses.
- so_what: prerequisite for *any* leakage claim; without it the paper's numbers are dismissible. A clean negative here is itself a contribution.

## T6 — Manifest grammar coverage / expressiveness (G1/G2)
- Hypothesis: ≤5 constructs cover >90% of benign pipeline reveals while vetoing all mutated adversarial variants — or fail on a nameable construct (e.g., adaptive burden-test definitions).
- MVE: encode artifact inventories of ~10 published multi-site GWAS pipelines as manifests; benign pass-rate + mutation veto-rate. Effort: weeks. Risk: LOW–MEDIUM.
- prior_work: declassification allowlists (Sabelfeld/Sands); policy-as-code; none with genomic reveal-object semantics.
- so_what: bounds how expressive pre-authorization must be — a design-defining diagnostic for the whole governance approach.

## T7 — Bypass-injection bake-off: enforced egress vs logging-only (proposal P0-5 / G4)
- Hypothesis: the MPC reveal gate blocks the bypass class by construction but has a measurable false-block rate; logging detects some with latency and silently misses others.
- MVE: instrument the same workflow three ways (SPU reveal gate / plaintext+Nextflow+in-toto logging / per-site API egress); inject identical bypass suite; measure caught/succeeded/latency/benign-blocked. Effort: months (scope to 2 instruments ⇒ 1–2). Risk: MEDIUM–HIGH.
- prior_work: SandboxEval escape taxonomy; provenance is observational; "MPC reveal = narrow enforcement point" never tested.
- so_what: directly falsifies/conserves the proposal's core backend hypothesis; a tie would honestly demote MPC to "one instantiation".

## T8 — Slicing attack: reconstructing forbidden matrices from authorized fine-grained reveals
- Hypothesis: scalar-level gates admit partial GRM/kinship reconstruction; vector/matrix gates don't; the threshold is predictable from rank algebra.
- MVE: reveal gate at scalar/vector/matrix granularity in SPU; attempt matrix reconstruction from authorized scalar streams under a query budget; reconstruction error vs granularity. Effort: weeks. Risk: MEDIUM.
- prior_work: repeated-query beacon attacks; DP composition (adaptive). No "granularity of the gate" attack on reveal governance itself.
- so_what: a named attack on the governance mechanism — the kind of adversarial result that makes a systems paper credible. Decides gate granularity.

## T9 — Statistical cost of pre-commitment (proposal P1-3 / G7)
- Hypothesis: discovery cost is small (few hits lost) while peek-then-adjust λ_GC inflation disappears; sign reverses for rare-variant analyses needing adaptive burden definitions.
- MVE: adaptive vs manifest-locked workflow on simulated panels with planted associations; compare power, top-hit overlap, λ_GC. Effort: weeks. Risk: MEDIUM.
- prior_work: pre-registration in clinical trials; none for multi-site GWAS statistical cost.
- so_what: quantifies whether governance destroys utility; a clean negative repositions the paper to "quantifying the cost of pre-commitment".

## T10 — Where composition additivity fails for genotype attacks
- Hypothesis: additive composition (PINQ-style reduction) holds for fixed panels but adaptive SNP selection concentrates risk super-additively at low-MAF SNPs.
- MVE: theorem for fixed-panel case; simulation adaptive vs fixed query selection; cumulative AUROC vs additive prediction. Effort: weeks–months. Risk: HIGH.
- prior_work: PINQ/Airavat/odometers (assume DP sensitivity); adaptive-query attacks.
- so_what: locates the exact boundary of what any composition rule could cover; thin theorem + honest violation map is a strong diagnostic.

## T11 — Sanitized-diagnostics paradox (inverted assumption)
- Hypothesis: *providing* sanitized diagnostics (binned spectrum, rounded λ_GC) reduces cumulative risk by removing analysts' incentive to probe — total-risk curves cross.
- MVE: scripted analyst policies doing QC under a fixed reveal budget, with/without sanitized diagnostics pre-provided; count reveals + cumulative attack success. Effort: weeks. Risk: HIGH (simulation-validity threat).
- prior_work: none — inverted framing of the standard "suppress" instinct.
- so_what: governance should sometimes *add* reveals; a counterintuitive, high-upside result if the simulation is credible.

## T12 — Receiver-dimension collapse: what GWAS artifacts actually become public?
- Hypothesis: most non-kinship artifacts become world-readable in supplements, collapsing the receiver dimension; expressiveness should go to object relations.
- MVE: systematic supplement audit of ~20–30 papers (pre-registered inclusion); tabulate artifact type × publication status. Effort: days–weeks. Risk: LOW.
- prior_work: none as a systematic audit.
- so_what: cheap, shapes the manifest grammar's design space; likely a citable negative on receiver-set governance.

## T13 — TEE-vs-MPC enforcement parity (executor lens: backend)
- Hypothesis: a TEE (simulated SGX) reveal gate enforces the same bypass class as MPC with a smaller TCB but a different trust assumption (hardware vendor); the *enforcement boundary* is comparable, the *trust story* is not.
- MVE: same bypass suite from T7 through an SGX-simulated gate vs SPU gate; compare caught/detected/out-of-TCB + TCB size. Effort: weeks. Risk: MEDIUM.
- prior_work: proposal only does a *qualitative* TEE comparison table (P2-3, "not implemented"); this makes it empirical.
- so_what: converts the proposal's weakest section (asserted TEE contrast) into measured evidence.

## T14 — Static reveal-compatibility decidability (executor lens: method-transfer)
- Hypothesis: for the restricted DSL, whether a (DAG, manifest) pair can ever produce a subset-consistency violation is decidable in linear time, and the decision procedure doubles as the compiler's static checker.
- MVE: encode the DSL/datagraph; implement the decision procedure; verify against P0-5's bypass suite (mutations should be rejected statically). Effort: weeks. Risk: MEDIUM.
- prior_work: IFC type systems (in-program); abstract interpretation; no workflow-level reveal-compat decision procedure.
- so_what: turns "unregistered reveal rejected at compile time" from a mechanism claim into a decidable property with a complexity bound — a small formal core for a systems paper.

## T15 — Enforce vs detect vs noisy-answer: head-to-head (executor lens: contradiction)
- Hypothesis: for honest-but-sloppy authors (A2), detection-only and enforcement produce identical leakage but different false-block rates; for the cross-run adversary (A3) only enforce+ledger helps. Locates where enforcement adds value over provenance.
- MVE: same bypass suite under detect-only / block / noisy-answer modes; (leakage, utility, false-block) frontier per adversary class. Effort: weeks. Risk: MEDIUM.
- prior_work: merges pi#7 logging comparison with the proposal's "detect don't block" cross-run choice; no three-mode head-to-head.
- so_what: answers "why not just log?" — the question a systems reviewer will always ask.

## T16 — Small-site safety boundary curve (executor lens: diagnostic)
- Hypothesis: there exists an exact (n_site, #SNPs, MAF floor, precision) boundary above which a released summary set is membership-inferable at TPR@FPR≥target, and it is reusable as a design curve.
- MVE: sweep n_site ∈ {50..2000} × SNP count × MAF × precision; measure TPR@FPR per cell; fit boundary. Effort: weeks. Risk: LOW–MEDIUM.
- prior_work: Homer/Sankararaman give single points; no operational boundary curve for multi-site site-size planning.
- so_what: directly actionable (tells sites when their stats are unsafe); distinct from T5's matched-control decomposition.

## T17 — Condition-bit leakage (executor lens: diagnostic; proposal P0-4 sub-item)
- Hypothesis: revealing *which branch* a secret-dependent decision took (the condition bit) leaks measurably more than the branch output itself for the PCA-count decision.
- MVE: same reveal set with/without logged condition bits; attack success delta. Effort: days–weeks. Risk: LOW–MEDIUM.
- prior_work: proposal's own definition 2 (branch decision is itself a reveal); never measured.
- so_what: tests whether the proposal's "condition bit counts as a reveal" rule is empirically justified or over-strict.

---

## Feasibility gate result
All T1–T17 pass (CPU-only, public data available, ≤ ~1 week compute). None dropped.
No dataset is provably unavailable; no idea exceeds budget.

## Dedup log
- pi#7 ∪ L4 → T7/T15 (kept both: T7 = instrument bake-off, T15 = adversary-class×mode head-to-head).
- pi#2 ⊃ L6-partial → T2 (condition-bit split out as T17).
- pi#5 ∪ L5 → kept separate (T5 decomposition vs T16 boundary curve).
