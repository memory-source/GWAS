# Senior ML Review — Loop 2, Round 2 (confirmation after fixes, 2026-10-04)

All round-1 findings re-verified directly from artifacts (changed code re-read in
full; both re-generated JSONs read verbatim; fresh pdflatex ×2; number-match of
the new s7 values; unchanged JSONs confirmed untouched by timestamp).

## Fix verification (round-1 → round-2)

1. Release-unit ambiguity → **FIXED in substance**: `release_objects` now
   returns `loadings`; s2c's T1/T3 consume `rel["loadings"]` (T2 correctly still
   fits the attacker's reference basis); re-run JSON is **bit-identical**
   (pcdist 0.8813, ΔT2 +0.0026±0.0047, ΔT3 +0.1020) — number unchanged as
   expected. Paper contribution 2 scopes the claim: "implementable from the
   release alone *if the study's PC loadings are shared*". RESIDUAL: the
   abstract still states "attacker-implementable from the release alone"
   unconditionally (one clause).
2. Operating point → **FIXED**: §5 Operating-point paragraph states
   900/900/900, 4000 sites (≈3814 QC), ratio ≈13, and seed counts (5 frontier/
   mechanism, 3 intermediate, 2 boundary/dim) — all verified correct against
   the JSONs' seed lists.
3. Uncertainty gaps → **FIXED**: s7 now averages 5 random three-way splits with
   stds; paper reports 0.684/0.681/0.564/0.640/0.808 (±0.009–0.014), control
   0.504, PC-distance 0.347 naive / ≈0.65 flipped — **all match the new
   s7_real.json**; ordering preserved (perpop > full/2dp > 1dp/top-k).
4. Text/numeric slips → **ALL FIXED**: "+0.003" (bogus ±0.006 gone); "0.84 at
   both 2dp and full precision" (0.8446 ✓); "same-signed but strongly
   attenuated" with both orientations and updated numbers; top-k criterion
   disclosed (|AF−refAF| proxy, no phenotype); Table 2 caption now cites §6
   for the 5-seed confirmation.
5. F/M taxonomy → **FIXED** (§6 defines Fn/Mn and "intermediate", explicitly
   handles per-site AF as intermediate-but-reported-separately). NIT: the Mn
   enumeration omits the condition bit (M4 in code; Table 2 row "+condition
   bit" carries no M-number).
6. Missingness caveat → **FIXED** (Discussion: structureless i.i.d. noise; real
   batch-correlated QC could differ).
7. Ratios rewording → **FIXED** ("largest at ratio ≈3–10, collapses toward
   extreme ratios where the AF attack alone is saturated") — matches s6 deltas.
8. `\label{sec:discussion}` → **FIXED**: fresh compile 0 errors, **0 undefined
   references**, 6 pages. 3 overfull hboxes remain (cosmetic, now disclosed).
9. Stale docs → **PARTIAL**: EXPERIMENT_RESULTS.md S7 section fully synced
   (new numbers + attenuated wording); tracker's "compiled clean" fixed; BUT
   the tracker's S7 table row still shows the old single-split key result
   (0.699/0.793/0.510 vs new 0.684/0.808/0.504). One line.
10. auroc() tie handling → unchanged, latent only (no reported number affected).

## Verdict

Score: **7.5/10** for a top venue (up from 7.0). Every substantive round-1
finding is fixed and verified; nothing measured is wrong; the remaining items
are one clause in the abstract, one tracker line, one optional word (M4), and
3 cosmetic overfull boxes.

**READY: Almost** — one clause from Yes. The two residuals: (a) add the
loadings conditional to the abstract sentence ("…from the release alone, if
the study's PC loadings are shared"); (b) update the tracker's S7 row to the
5-split numbers. Both are one-line edits; no re-experimentation remains.