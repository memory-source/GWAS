# Cross-model review verdict — round 4, final revision check (2026-10-02)

Reviewer scope: PoPETs/USENIX. I verified each checklist item against
`paper/main.tex`, the raw JSONs, and the code (mtimes + content). I also
externally verified the bibliography entries that could not be confirmed from
the repo. This is the final verdict.

## 1. Checklist verification (all nine items)

| # | item | verified? | evidence |
|---|---|---|---|
| 1 | Abstract reworded (intermediate vs single channel) | **YES** | Abstract now: PC diagnostic = "largest intermediate-release leakage channel"; per-site AF = "largest single channel of all (0.86)". The (ii)/(iv) contradiction is resolved and both statements match the data (0.791 PC-alone, 0.858 per-site, 0.758 base) |
| 2 | msprime cross-check reported | **YES** | Boundary section reports 0.642 full / 0.649 1dp, "quantization effect is within noise" — numbers match `s4_msprime.json` exactly (0.6420/0.6488); framing is honest and consistent with the saturation reading |
| 3 | Table 2 uses 5-seed S2c values | **YES, with one nit** | 0.758/0.860/+0.102 match `s2c_mechanism.json` (0.7576/0.8598/0.1022). Nit: the +λGC/+condition-bit/+missingness rows derive from the unchanged 3-seed `s2b` run (mtimes 04:06/04:07); the Δ values are seed-robust (0/0/+0.001), but the caption says "5 seeds" without provenance, and the +missingness row is internally inconsistent (AUROC 0.758 with Δ+0.001 ⇒ 0.759) |
| 4 | T3 range fixed | **YES** | "+0.08 to +0.11, mean +0.10" matches per-seed T3 deltas (0.083–0.112, mean 0.1025) |
| 5 | Permutation control added to s2c | **YES** | `s2c_mechanism.py` lines 90–109 (T5, within-ancestry permutation, added 04:54); JSON: per-seed [0.520, 0.506, 0.491, 0.463, 0.508], mean 0.498, std 0.020 — paper's "0.498±0.020" matches exactly |
| 6 | AF release relabeled | **YES** | "0.75 at 2dp; 0.72 at full precision" in abstract, contribution 1, and the mechanism paragraph — consistent with s1c (F3 0.714/F5 0.720) and s2c base 0.758 at 2dp |
| 7 | Dimension-sweep sentence reworded | **YES** | "nine operating points (ratios 1.7–60), always positive (+0.03 to +0.12) ... beats the AF summary release at all nine points" — matches `s6_dim.json` exactly. Minor: "largest at small sites" — the single max (0.124) is at n=100 but the mean Δ is largest at n=300 (0.103 vs 0.076/0.078); "small-to-moderate sites" would be tighter (optional polish) |
| 8 | Real bibliography (15 entries) + concrete A1/A2/A3 | **YES** | 15 entries; externally verified: homer2008, sankararaman2009, shringarpure2015, bu2021, ostrak2021, wagner2016, intoto, rfc6962, prov-o are standard; **ppgwas2025 verified real** (Nat Comms, s41467-025-66771-z), **igwas2023 verified real** (PoPETs 2023-0026), **fedlmm2025 verified real** (bioRxiv preprint), **popets2026 verified real** (PoPETs 2026(4):542–560, popets-2026-0134, Berrang/Ryan/Wooldridge — and its content matches how the related-work section uses it). A1/A2/A3 each now have a concrete capability sentence. Three label/hygiene nits remain (see §3) |
| 9 | Self-reference demoted; contributions lead with the finding | **YES** | Contribution 1 is now the diagnostics finding; the correction is a "Methodological note" paragraph in the protocol section, framed as "one of the paper's methodological takeaways" |

Every number in the paper was re-checked against the raw JSONs: Table 1 = s1c
(including the per-ancestor caption), Table 2 = s2c, boundary table = s5,
mechanism paragraph = s2c (pcdist 0.791, pcalign 0.504, T2 +0.003±0.006, T3
+0.1025±0.0104, eig ratio 0.9998), discussion S3b = s3b (3→0 both families,
equal measured quality, small AUROC change). No unsupported numbers remain.

---

## 2. Final verdict

**YES — submission-plausible as a measurement paper.** All nine round-3
checklist items were genuinely applied (verified in the paper, the code, the
data, and externally for the citations). The paper is internally consistent,
every number is backed by an artifact, the mechanism is identified with proper
controls (independent-spectrum collapse + attacker-realistic centering intact),
the negative results are reported without spin (λ_GC and condition bit at
exactly zero; S3b's measured AUROC change described as small), the limitations
are plain, and the stated next step (real released AF vector, e.g. 1000
Genomes populations as pseudo-sites) is the correct one. What remains is
text-level hygiene only (§3).

---

## 3. Residual must-fix items (text only)

1. **popets2026 bibitem label** says "Wooldridge et al.(2026)" — the first
   author is Berrang (P. Berrang, M. Ryan, K. Wooldridge). Relabel to
   "Berrang et al.(2026)".
2. **fedlmm2025 entry has no author names** (label says "Chunduru et al.(2025)"
   but the entry body lists none). Add the real author list from the bioRxiv
   page (2025.12.16.693409).
3. **igwas2023 bibitem label** says "Pascoal(2022)" while the entry cites
   PoPETs 2023(2) — align the label year to 2023.
4. **Table 2 caption provenance**: λ_GC/condition-bit/missingness deltas derive
   from the 3-seed s2b run (seed-robust: 0/0/+0.001), while final-only and +PC
   are 5-seed S2c. Add a one-line caption footnote, and fix the +missingness
   row's rounding (0.758 with Δ+0.001 ⇒ 0.759).

Optional polish (not blocking): "largest at small sites" → "largest at
small-to-moderate sites"; the natbib `[numbers]` option mixed with author-year
`\bibitem` labels is an unusual combination — pick one citation style for
camera-ready.

---

## 4. One-line summary for a PC member

A rigorous, fully simulation-based measurement paper showing that in multi-site
GWAS the released study-fitted PC spectrum — not the summary statistic — is the
dominant intermediate membership-leakage channel (0.79 alone, 0.86 combined,
vs 0.72–0.76 for AF releases), with the mechanism identified as study-fit
in-sample-ness under clean disjoint-negative controls and honest, plainly
stated limits; novel question and reusable methodology, with real-data
validation as the named next step.
