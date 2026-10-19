# Round-5 Re-Review (autonomous review loop round 3/3) — 2026-10-04

Reviewer: same independent artifact-level verification. Read the revised main.tex in
full, read code/s15_controls.py and results/raw/s15_controls.json, re-ran the full
s15 experiment from seeds (bit-identical, diff 0.00e+00 on all 9 batch numbers and all
4 structure numbers), recompiled (8 pp, 0 errors, 0 undefined refs), truncation sweep
clean.

## Fix verification (all four confirmed)

1. **Line 254 sign error FIXED**: "The train-calibrated orientation is the
   members-farther direction on 25/25 seeds" — matches s10_paired.json (pc_sign=+1.0,
   sign_agreement 1.0) and the Mechanism paragraph. Residual nit: "agreeing with the
   mechanism paragraph below" — the Mechanism paragraph is ABOVE (line 233 vs 248);
   "below" should be "above". Cosmetic direction word, not a data error.
2. **"large-panel" regression FIXED**: "sparse-marker/small-panel corner (400
   individuals, 500 markers)" — consistent with the S12 grid (400 = smallest panel,
   500 = fewest markers) and the data (peak +0.163 at n400_m500).
3. **Batch gloss FIXED**: "(0.76, below the 0.85 full-coverage AF attack and the 0.80
   uniform-random-70%-coverage control on the same panels; s15_controls.py)" — all
   three numbers now artifact-backed (s15: 0.7592 / 0.8458 / 0.8004), bit-identical
   on re-run.
4. **Fst estimator named and backed**: paper cites "≈0.10 vs ≈0.05" (s15:
   0.0994/0.0503) and "0.26 vs 0.14" (0.2583/0.1422) with "both computed in
   s15_controls.py". The estimator is named in the s15 docstring (Var_k(p_k)/(p̄(1−p̄)),
   ddof=0) and the code's Nei G_ST form (Ht−Hs)/Ht is algebraically identical (I
   verified: Ht−Hs = 2·Var_ddof0(p_k)). Direction (real MORE structured) robust under
   ddof=0 and ddof=1.

## Score: **8/10** (up from 7.5)

Every cited number now traces to results/raw/*.json with bit-identical reproducibility
(S9c and S15 re-runs: 0.00e+00 diff). Claims correctly scoped; the previously false
explanation replaced by an artifact-backed, correctly-signed statement plus a flagged
hypothesis; the real-data increment has a formal paired test; the batch caveat is
measured and correctly interpreted; the batch gloss and Fst numbers are artifact-backed.
Caps below 9 are the disclosed residuals: one real dataset (chr22), headline magnitudes
simulation-only, PC channel tested with one attack family, no msprime cross-check for
the PC channel, no MPC prototype — all honestly listed as limitations/future work.

## READY? **Yes**

Submission-plausible as a measurement paper. One cosmetic word ("below" → "above" in
the paired-significance paragraph) recommended before submission; it references no data
and cannot mislead about any result.
