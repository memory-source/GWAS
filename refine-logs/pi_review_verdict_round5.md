# Cross-model review verdict — round 5: real-data design + second attack family (2026-10-03)

Reviewer scope: cross-family review of two approaches before finalization. I
read `code/s7_real.py`, `code/s8_two_family.py`, both s8 JSONs, the rewritten
`common.py` (mtime Oct 2 17:39), and verified the arithmetic empirically.
**A third, critical finding emerged (§3.0) that outranks both questions: the
maf() rewrite breaks the reproducibility of every number in the paper.**

---

## 0. Critical finding first: the maf() rewrite invalidates all paper tables

`common.py`'s `maf()` was rewritten Oct 2 17:39 — **after** every run that
produced the paper's numbers (s0/s1c/s2c/s3b/s4/s5/s6, all 04:xx) and **before**
s8 (Oct 3 00:49).

- Old formula `min(G.mean(0), 1 − G.mean(0))` goes **negative** for any SNP with
  alt-count frequency > 1 (AF > 0.5), so the old qc_mask [0.05, 0.95] silently
  **dropped every SNP with AF > 0.475**. Verified on seed 0 data: old mask kept
  **1866/4000 SNPs**, all with AF ≤ 0.475; **2110 common SNPs (52.75%) were
  excluded** from every "all tested SNPs" tier.
- The fixed formula (`min(af, 1−af)`, af = alt-count/2) keeps **3812/4000**
  (AF ∈ [0.05, 0.95]) — including all common variants.

Consequences:
1. **The current pipeline does not reproduce any number in the paper.** s8, run
   with the fixed maf and a nominally identical config (n=1000/group,
   4000 SNPs, 300/300/300 split, seeds 0–4), gives AF-LRT 0.840 at 2 dp and
   0.845 at full precision — versus the paper's Table 2 base 0.758 and Table 1
   F5 0.720. Same seeds, same nominal config, different data: the release
   definition changed.
2. **The paper's release tiers were measured on a biased release.** Tables 1–3
   and the mechanism/boundary/sweep sections describe tiers like "all SNPs, AF
   1 dp / 2 dp / full" — but the actual tested set excluded every common
   variant. Common SNPs carry substantial AF-LRT signal, so the corrected
   numbers are **higher**: the paper's claims are directionally conservative
   (real leakage is worse than reported), but every table is stale.
3. The bundle's "run on the SAME panels as the AF-LRT" claim for s8 is false
   relative to the paper: within s8 all families share panels, but s8's panels
   were built with the fixed maf and are not comparable to Table 1/2 as-is.
4. The tier ordering and every qualitative conclusion survive in s8
   (F4-1dp 0.62 ≪ F3-2dp 0.84 ≈ F5 0.85; PC-distance 0.88 tops all three
   families at every tier, every PC seed above every AF-LRT seed) — the story
   is intact; the numbers must be regenerated.

**Required action (not optional):** re-run s0–s6 with the current fixed maf,
regenerate every table and figure, update all in-text numbers, and state the
MAF filter actually applied in the method section. Pin the pipeline version in
the repo (the paper and code must correspond). Note the direction honestly:
the corrected filter makes the leak larger, which strengthens the paper.

---

## 1. (A) Are F2/F3 legitimate distinct MI families, or re-parameterizations of F1?

**They are re-parameterizations of the same underlying signal — say so plainly.**
All three are functionals of one quantity: individual-vs-AF proximity.
- F1 (AF-LRT): linear in the genotype, weights = per-SNP log-likelihood ratio.
- F2 (AF-corr, verified implementation: per-row Pearson correlation with the
  study-AF vector minus correlation with the ref-AF vector): linear in the
  genotype, weights = the centered AF vectors.
- F3 (ML-dev, verified: CV logistic on two features — mean squared deviation
  from study AF and from ref AF): quadratic functional of the same proximity,
  with a learned monotone combination.

None exploits information the others do not (no haplotype/LD-phase structure,
no rare-variant burden, no query-pattern structure). The s8 docstring's "three
structurally distinct membership-inference families" overclaims, and the paper
must not repeat it. The correct claim: **"three functional forms of the
AF-proximity attack agree on the frontier ordering, and the PC-diagnostic
channel tops all of them"** — which the data does support (PC 0.882 > all
three at every tier).

What the exercise genuinely buys: (i) robustness — the result is not an
artifact of the LRT's specific functional form; (ii) one real minor finding:
**at 1 dp the correlation/deviation forms beat the LRT** (0.660/0.659 vs
0.622) — the LRT is more fragile to coarse rounding than correlation-style
scores; worth one sentence in the boundary section, since it says the
coarse/fine lever interacts with the attack's functional form. Removing the
beacon was right — it was null by construction, and `s8_two_family.json`
confirms it (0.489 ≈ chance) before its replacement.

If the "one attack family" clause is to be closed *properly*, the missing
family is a channel that reads different information — a haplotype/LD-based
attack (Bu et al. 2021, already in the bibliography) is the natural one. That
is optional work, not a requirement of this review; the honest reframing is
the required part.

Provenance nit: `s8_families.json` has no config block (no n_site/n_snps/seed
record), and the guard comment in `s8_two_family.py` says the ML features are
tier-independent while the code passes the tier's AF into `ml_dev_features`.
Add a config dict; fix or implement the comment.

---

## 2. (B) Correctness risks in the S7 design

**S7 has not produced results** (no `s7_real.json`; the 76 MB chr22 VCF and the
panel are downloaded, the script written Oct 3 00:45). Reviewing the approach:

1. **FATAL — VCF GT parsing bug.** `p = line.split("\t", 9)` makes `p[9]` the
   *entire sample block as one tab-joined string* (maxsplit=9). `gt[i]` then
   character-indexes that string — `idx` holds *column positions* (~9…2512), so
   `gt[i]` reads an arbitrary character (an allele, a separator, or a tab of
   some other sample), not sample i's genotype. Also `b = g[2] if len(g) > 2`
   never fires (`g` is a single character), so "genotypes" are binary {0,1},
   never 2. Every downstream quantity would be an artifact of character
   positions. Fix: split the sample block fully (`p[9].split("\t")`), index by
   the column offset (idx − 9), parse both alleles of the 3-char GT field, and
   assert no "." alleles (1000G Phase 3 GTs are complete — assert it).
2. **Split-ordering confound.** `three_way_split` carves contiguous
   panel-file-order thirds within each super-population. 1000G panel order
   groups subpopulations (CEU/GBR/TSV/… under EUR), so the study, reference,
   and non-member panels can differ in subpopulation composition — an ancestry
   confound *within* super-population that the ancestry control would then
   correctly flag but the design should not have. The simulator's contiguous
   slices were fine (i.i.d. individuals); real data is not exchangeable in file
   order. Fix: shuffle indices with a recorded seed (or stratify the split by
   subpopulation) before carving thirds, and run a few random re-splits to
   report panel-assignment variability (the current bootstrap covers individual
   resampling only).
3. **"LD-thin" is distance-thinning.** `ld_thin` keeps one SNP per ≥100 kb
   window — a physical-distance proxy, not r² pruning. Kept SNPs can still
   share long-range LD, and the summed-LR attack assumes SNP independence, so
   retained LD inflates it. Either rename ("distance-thinned, ≥100 kb") or do
   proper r² pruning (bcftools/plink) — rename is the minimum.
4. **F2_topk proxy uses the attacker's reference.** Top-200 by |AF_s − AF_r|
   defines the "released interesting set" from the attacker's own panel; that
   conflates release policy with attacker knowledge. Name it explicitly as a
   proxy (no phenotype exists in this data), or select by a release-side
   criterion.
5. **Correct as designed (verified):** per-population AF release (F5_perpop)
   matches the paper's per-site framing; the two controls are the right pair —
   `structure_only_control` on a *reference-fitted* basis (both target panels
   out-of-sample, expected ≈0.5) versus `study_pc_distance_auroc` on the
   *study-fitted* basis (the S2c leak, expected >0.5); no phenotype needed for
   the MI tiers. The new 4-arg `anc_only_attack` signature matches the
   rewritten `common.py`.

---

## 3. Does (A)+(B) close the rejection clauses?

**(A) — partially, and only with a reframe.** What closes: "the result is not
an artifact of one attack's functional form" (three forms agree; PC tops all).
What does not close: "structurally distinct attack families" — F2/F3 are
re-parameterizations of AF proximity. The clause becomes "one attack *channel*
(AF proximity, three functional forms) plus the PC-diagnostic channel." A
haplotype-based attack would close it fully; the reframing is mandatory either
way.

**(B) — not yet.** The design is right and the data is downloaded, but the
script has a fatal parsing bug and the split-ordering confound, and it has not
run. The "all simulation" clause closes only when S7 runs correctly and
reports honestly — chr22, three pseudo-sites, no phenotype is a limited but
legitimate external-validity check if scoped as such.

**Still missing, concretely (no invented requirements):**
1. Regenerate every paper number with the fixed maf (§0) and make paper and
   pipeline correspond; disclose the filter actually applied. This is the
   blocking item — nothing else should be finalized on stale tables.
2. Fix S7 (parsing, randomized/stratified split, honest thinning name, GT
   assert), run it, and report it in the paper with its stated scope.
3. Reframe the s8 claim (functional forms, not distinct families); add a
   config block to `s8_families.json`; fix or implement the guard comment.
4. Re-run the ancestry-only control with the new 4-arg function: the old 2-arg
   path (as called in s1c) never scored study members — the new docstring says
   so itself — so Table 1's caption "ancestry-only control at chance at every
   tier" is currently a panel-identity statement, not a member-vs-non-member
   ancestry attack. With the new function on ancestry-matched panels it should
   stay ≈0.5; make the caption true as written.

---

## One-line summary

The approaches are directionally right — three functional forms of the
AF-proximity attack agree and PC still tops them, and the 1000G design is
structurally sound — but the maf() rewrite means every published number is
stale (the old pipeline silently dropped all 2110 SNPs with AF > 0.475), S7
cannot run as written (GT parsing bug + split-ordering confound), and the
"distinct families" claim must be downgraded to "functional forms of one
channel" before anything ships.
