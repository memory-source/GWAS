# RESULT-TO-CLAIM EVALUATION (Step 2)

**Claim under test:** PC diagnostic is the dominant intermediate membership-leakage channel (~0.88 > AF summary ~0.84, +0.10 increment); per-site AF release is the largest single channel (~0.96); mechanism = in-sample-ness of the study-fitted spectrum; disjoint-negative shows leakage is real; ordering reproduces on real 1000G chr22; self-reference error inflates ~0.13.

**Verification:** all 14 cited values re-derived from `results/raw/*.json` (not taken from notes). Internal consistency holds: S1d `disjoint_mean` 0.8446 == S1c F5_full exactly; S2c base 0.8392 ≈ F3(2dp) 0.8396.

## 1. claim_supported: **partial**

## 2. What the results support
- **PC > AF summary in simulation — and the evidence is paired, not cross-script:** S8 same-seed/same-panel comparison at the F5_full tier: PC 0.8816 vs AF-LRT 0.8446, **5/5 seeds** PC wins (paired diffs +0.023…+0.049, mean +0.037). S2c pcdist_only 0.8813±0.0053.
- **+0.10 increment:** S2c base+M2 delta +0.1019±0.0074; combined 0.9411 is +0.0966 over the full-precision AF release 0.8446. Both round to ~+0.10. (S2b independently: M2 marginal +0.098; M1_λGC and M4_condition-bit marginals exactly 0.)
- **Per-site largest single channel:** 0.9642±0.0037, per-ancestor 0.961–0.969, TPR@1% 0.588. Claim correct.
- **Mechanism = in-sample-ness (best-supported sub-claim):** T2 (spectrum refit out-of-sample) kills the increment: +0.0026±0.0047, sign flips across seeds → statistically ≈0. T3 (released spectrum + attacker-side study-mean centering) preserves it: +0.1020±0.0077. pcalign-only 0.5012 (chance → not ancestry alignment). Permuted labels 0.486. Released/in-sample eigenvalue ratio 0.9998. Five convergent controls.
- **Not an ancestry confound:** ancestry-only 0.5019 (S1c), structure-only 0.5036 (S7 real), pcalign 0.5012, permuted labels 0.486, disjoint non-member panels throughout.
- **Self-reference inflation ~0.13:** 0.9765 vs 0.8446 = +0.1319 (per-seed 0.119–0.143). Verified.
- **Boundary S5:** n=50/site, 2dp → AUROC 0.9995, TPR@1% 0.9933 (rounding does not protect at small n); at n=800 the signal is weak regardless.
- **"Families are one channel":** S8 at F5_full: f1 0.845 / f2 0.833 / f3 0.832 — within 0.013 of each other, all below PC 0.882. Consistent with functional forms of one information channel.
- **AF-tier ordering reproduces on real chr22:** perpop 0.8080 > F5 0.6839 > control 0.5036 (5 splits, ±0.012–0.013).

## 3. What the results do NOT support
- **"The ordering reproduces on real chr22" as stated.** The PC-dominance part does not reproduce: real PC-distance raw 0.3472 (anti-predictive); even generously sign-flipped ≈0.653 **<** full-AF 0.684. Only the AF-release tier ordering (per-site > full > control) reproduces. The paper's own S7 artifact contradicts the full-ordering reading.
- **Calibration-dependence of the 0.88:** at fixed orientation the raw diagnostic is anti-predictive even in simulation (S1c `study_pc_distance_auroc` 0.1184). The 0.8813 requires attack-side sign calibration (S2c CV meta-attack learns it). Legitimate for an adaptive attacker, but the claim "projection-distance attack reaches 0.88" must state this operating detail.
- **Operating-point generality of "+0.10":** S6 sweep shows pcdist 0.69–0.99 and delta 0.006–0.136 across the (n_snps × n_site) grid. The headline values are specific to ratio ≈ 13 (300/site, 4000 SNPs); at ratio 60 (6000 SNPs, 100/site) the increment is +0.006. The claim as worded reads as general.
- **"Dominant intermediate":** "intermediate" is undefined; the PC-vs-AF gap (+0.037 paired) is real but modest — "dominant" overstates 0.88 vs 0.845 unless scoped.

## 4. Missing evidence
1. Sign-calibrated / oracle-orientation real-data PC attack + simulation at the 1000G-matched operating point (450/panel, 1363 SNPs, 15 kb thinning) to localize the attenuation (n? SNP count? LD/thinning?).
2. Formal paired significance test for PC vs AF (5/5 seeds is strong but unquantified — paired t/Wilcoxon or paired CI).
3. Cross-dataset replication: currently one chromosome, one VCF, 5 splits.
4. Explicit release-tier taxonomy defining "intermediate" in the claim text.
5. Batch-correlated missingness variant for M3 (current missingness structureless; marginal +0.0004 may not hold under realistic batch structure).
6. Already-flagged, correctly labeled: no MPC/governance prototype (argued, not built); S3b is a scripted-analyst conjecture.

## 5. Suggested claim revision
> "At the paper's operating point (3 ancestries × 300 individuals/site, 4000 SNPs), the released study-fitted PC diagnostic is the largest non-allele-frequency membership channel: a projection-distance attack with attack-side sign calibration reaches 0.881±0.005 AUROC, above the full-precision allele-frequency summary release (0.845±0.009; paired per-seed comparison 5/5 seeds, mean gap +0.037), and adds +0.102±0.007 on top of it. The per-site (per-stratum) allele-frequency release is the largest single channel (0.964±0.004). The mechanism is the in-sample-ness of the study-fitted spectrum: refitting the spectrum out-of-sample removes the increment (+0.003±0.005), while replicating study-mean centering on the released spectrum preserves it (+0.102±0.008); released eigenvalues match in-sample ones (ratio 0.9998) and alignment-only and permuted-label controls sit at chance. Under disjoint-negative evaluation the individual-level leakage is real (structure-only controls 0.502 sim / 0.504 real). A self-reference control-design error inflates naive estimates by 0.132 AUROC (0.976 vs 0.845). On real 1000 Genomes chr22 the allele-frequency tier ordering reproduces at smaller magnitude (per-population 0.808 > full 0.684 > control 0.504); the PC channel does not cleanly reproduce (naive 0.347, sign-flipped ≈0.65, below the full-AF release), so PC dominance is supported in simulation only."

## 6. Next experiments needed
1. **Real-data PC orientation study** (naive / sign-flipped / oracle) **+ matched-operating-point simulation** — the single most important experiment; explains the sim→real gap or bounds PC dominance to simulation.
2. **PC-vs-AF crossover map:** vary n_per_panel at matched SNP count/thinning; find where PC exceeds AF (if anywhere) → turns "dominant" into a scoped, defensible statement.
3. Second real dataset (chr1, or HGDP) replication of the AF-tier ordering.
4. Batch-correlated missingness variant of M3.
5. Optional future work: haplotype/LD-based attack as a genuinely distinct channel.

## 7. confidence: **high**

Numbers are deterministic artifacts, all verified; internal cross-checks (S1d==S1c disjoint, S2c base≈F3-2dp) pass; the failing sub-claim (real-data PC ordering) is directly contradicted by the paper's own artifact, and the user's caveats already acknowledge it.