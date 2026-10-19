# Round-5 fresh review — Round 3/4

Reviewer: pi (same method). main.tex re-read in full; make_figures.py fig7a and
s26_ibd.py docstring checked; residual-claim grep; fresh compile (12 pp, 0
errors, 0 undefined refs).

## Verification of the four claimed fixes — all applied

1. **Crossover restatement: APPLIED.**
   (a) §7 line-544 sentence rewritten to the measured-terms statement ("crosses
   zero at a measured structure statistic ≈0.06 — below the real 1000G level, so
   the real panels sit beyond the crossover"); the old nominal-0.10 claim is
   gone; the paragraph is no longer self-contradictory.
   (b) Abstract and contribution 5 now read "crosses at measured structure
   ≈0.06, below the real-data level" — consistent with §7 and with the s21 JSON
   (measured crossover ≈0.062 at n800; ≈0.043 at n1600; real level ≈0.10).
   (c) Fig 7a caption corrected the same way.
   (d) The shade was actually moved: axvspan(0.06, 0.10), label "crossover →
   real 1000G", and fig7_governance.pdf regenerated after the code change
   (mtime 03:17 > 03:16).
2. **"Phased" residue: APPLIED.** §7 heading → "An IBS window-sharing channel";
   Fig 7d caption → "An IBS window-sharing channel"; Limitations → "an IBS
   window-sharing proxy channel". Only line 480's "full phased-haplotype/IBD
   attack (Bu et al.) remains future work" retains "phased" — correct usage.
3. **Discussion sentence: APPLIED.** "which our ledger records (from
   producer-asserted metadata) as a refinement relation" — consistent with §7.
4. **s26 docstring: APPLIED.** "an 8-marker WINDOW-SHARING count … (unphased
   dosages)"; "real_phased.py loads the 1000G haplotypes but this arm does not
   yet exploit H0/H1 separately". The "longest shared run" and "real PHASED
   1000G haplotypes" claims are gone.

Residual-claim grep is clean: no "crosses zero near", no "structure level of
real 1000G panels", no "needs/requires discrete", no "A phased/IBS", no
"phased/IBD proxy", no "longest shared run". Regression check: headline numbers,
bibliography, governance numbers, canonical scoping — all unchanged and hold.
Compile clean (12 pp).

## Score: 9/10

## Verdict: Almost — one figure-axis detail remains

## Remaining weaknesses + minimum fixes

1. **Fig 7a band vs its own axis (the last item, one line).** The band now sits
   at axvspan(0.06, 0.10) on the NOMINAL-Fst axis, and the caption says it
   "marks 0.06–0.10, i.e. the crossover and the (higher) real-1000G structure
   level". Those endpoints (0.06 measured crossover, 0.10 measured real level)
   correspond on the nominal axis to ≈0.09 and ≈0.16 (measured ≈ 0.63 × nominal,
   from s21's fst_like). So the band as drawn covers the nominal crossover
   region (0.067–0.09) but NOT the real level (nominal ≈0.155) — the "→ real
   1000G" label and the caption's "i.e." clause overstate what the band
   contains. This is the same unit-mismatch class, now confined to the figure.
   *Minimum fix (one line + caption tweak):* draw the band at nominal
   ≈0.09–0.16 (the nominal equivalents of measured 0.06–0.10) with the caption
   saying "shown at their nominal-Fst equivalents", or add a secondary
   measured-structure axis.

2. **Cosmetic (non-blocking).** (a) §7 states the measured crossover as "≈0.06"
   (the n800 value); n1600's is ≈0.043 — "≈0.04–0.06" would be tighter (the
   direction, below the real level, holds either way). (b) The rewritten first
   sentence of §7 and the "Note the crossover is best read in measured structure
   terms" sentence are now redundant — both state the same measured crossover;
   merging them would tighten the paragraph. (c) "decreases monotonically at
   every site size" in the headline sentence vs the in-paragraph qualification
   ("strict at n=800 and 1600; at n=300 near ceiling, within tolerance of flat")
   — acceptable as-is since the qualification is present.

## Assessment

All substantive claims in the text are now correct and internally consistent:
the s22 interpretation matches the experiment's parameterization, the crossover
is stated in measured terms everywhere, the IBS channel is described as what
the code computes, the ledger's flags are honestly labeled producer-asserted,
and the paradox paragraph carries both the corollary and the negative-side
control. Every checked number matches the raw JSONs. The only remaining defect
is the Fig 7a band's axis units — a one-line fix. After it, submission-ready.
