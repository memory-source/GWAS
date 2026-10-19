# Round-5 fresh review — Round 2/4

Reviewer: pi (same method as round 1). main.tex re-read in full; s26 docstring,
make_figures.py, tracker checked; every changed claim verified against
results/raw/*.json; fresh compile (12 pp, 0 errors, 0 undefined refs).

## Verification of the six claimed fixes

1. **s22 inversion (W1): FIXED.** Verified: no "requires/needs discrete
   structure" language anywhere (only corrective occurrences). §7 now states the
   sweep "does not isolate discrete-vs-continuous structure" and "tracks
   structure strength rather than the discrete/admixed label" — accurate.
   Numbers match s22 JSON (−0.1147→"−0.12" most concentrated; +0.0957→"+0.10"
   most uniform; "nearly-pure individuals and near-uniform mixes" is the correct
   Dirichlet reading). Abstract and Limitations updated ("monotone in measured
   structure"; "beyond the crossover, not anomalous").
2. **Crossover conflation (W2): PARTIALLY FIXED — regression in summary
   positions.** §7's corrective note is present and correct (measured crossover
   ≈0.06, below real ≈0.10, real panels beyond; strict at n800/1600; n300 near
   ceiling, flat within tolerance). BUT the ORIGINAL conflated sentence at
   line 544 was not struck ("crosses zero near F_ST≈0.10 … the structure level
   of the real 1000G panels") and now directly contradicts the corrective note
   three sentences later in the same paragraph. The abstract (line 65) and
   contribution 5 (line 135) were reworded to "crosses at the real-data level" —
   which re-asserts the conflation at the two most-read positions and contradicts
   the corrected §7. Fig 7a caption (line 600) still says "crosses zero at the
   structure level of real 1000G panels (shaded)", and the Fig 7a shade was NOT
   repositioned (axvspan(0.09, 0.15) unchanged; fig7 PDF regenerated at the old
   setting) despite the change list claiming it.
3. **IBD not phased (W3): CORE FIXED.** §7 body reworded correctly ("IBS
   window-sharing … run here on unphased dosages … real_phased.py loads the
   haplotypes, but this arm has not yet exploited H0/H1 separately"); abstract
   and contribution 5 now say "IBS". Residual "phased" labels at three secondary
   spots: §7 heading "A phased/IBS channel" (582), Fig 7d caption "A phased/IBD
   channel" (603), Limitations "a phased/IBD proxy channel" (637) — inconsistent
   with the body's honest disclosure. (Line 480's "full phased-haplotype/IBD
   attack remains future work" is correct usage.)
4. **Anti-composition (W4): FIXED in §7.** "producer-asserted metadata (the
   prototype records and audits relations, it does not content-infer them)" ✓.
   Minor: the Discussion still says "our ledger's relation checks flag as a
   refinement" without the qualifier.
5. **s25 corollary + null side (W5): FIXED.** Δ_E=e(p_stop)−e(p_san), families
   set only the stopping rung, "parity is a consequence, not a coincidence", and
   the null family at p_san≥2 with negative Δ_E≈−0.24 ("a coarser-than-needed
   sanitized block IS the leak") — all present; −0.2437/−0.2530 match the JSON.
6. **Minor (W6): FIXED except the shade.** Seed inventory now covers the Fst
   sweep (15/cell) and admixture (10/setting); "loadings≫eigenvalues at every
   point" dropped and re-pointed to the real-data granularity section (the
   retained "own-basis ≫ other-basis at every point" IS measured per cell);
   s26 docstring's 15kb-thinning control line removed; tracker fixed to
   12 pages. Residual: s26 docstring still describes "longest shared run of
   identity-by-state" and "real PHASED 1000G haplotypes" (implementation is an
   8-marker window count; the real arm is unphased).

Regression check: headline numbers (0.88/0.84/0.94/25-25, 0.96, 0.99/1.00,
+0.13), bibliography (Berrang 2026, Suo 2025), canonical scoping, msprime
honesty — all hold. Compile clean (12 pp).

## Score: 8.5/10

## Verdict: Almost

One must-fix remains (text-only, 4 spots + 1 figure detail): the crossover
conflation survives at summary positions and now self-contradicts the corrected
§7 paragraph. Everything else verified correct; no data impact anywhere.

## Remaining weaknesses (ranked) + minimum fixes

1. **Crossover restatement regression (must fix).** (a) Strike or rewrite the
   line-544 sentence to the measured-terms statement; (b) abstract line 65 and
   contribution 5 line 135: "crosses at the real-data level" → "crosses at
   measured structure ≈0.06, below the real ≈0.10 — the real panels sit beyond
   the crossover"; (c) Fig 7a caption line 600: same correction; (d) reposition
   or relabel the Fig 7a shade (e.g., "measured 0.06–0.10: crossover and real
   level") to match the axis it sits on, or drop the shade claim from the
   caption.
2. **"Phased" label residue (minor).** §7 heading (582), Fig 7d caption (603),
   Limitations (637): "phased/IBS" → "IBS window-sharing"; "phased/IBD proxy
   channel" → "IBS proxy channel".
3. **Discussion composition sentence (minor).** "relation checks flag as a
   refinement" → "record (producer-asserted) a refinement relation", consistent
   with §7.
4. **s26 docstring (minor).** "longest shared run of IBS" → "8-marker
   window-sharing count"; "real PHASED 1000G haplotypes" → note the real arm
   runs on unphased dosages (H0+H1 collapsed).

## READY? Almost — one consistent-wording pass away.

No re-runs needed; all fixes are text-level plus one figure label. After item 1
(ideally 2–4), submission-ready.
