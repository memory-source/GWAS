# Cross-model review bundle — round 4 (revision check)

You reviewed round 3 and returned "submission-plausible, conditionally" with a
9-item checklist. The executor applied all items. Verify and give a final
verdict.

## Applied checklist
1. Abstract reworded: PC diagnostic = largest *intermediate* channel; per-site AF
   = largest single channel (fixes the (ii)/(iv) contradiction).
2. msprime cross-check now reported in the boundary section: 0.642 full / 0.649
   1dp, framed as "at that operating point the quantization effect is within
   noise," consistent with the saturation reading.
3. Table 2 now uses the 5-seed S2c values (base 0.758, +PC 0.860, Δ+0.102).
4. T3 range fixed to "+0.08 to +0.11, mean +0.10".
5. Label-permutation control added to `s2c_mechanism.py` and reported
   (0.498±0.020).
6. "full AF release (0.75)" relabeled to "0.75 at 2dp; 0.72 at full precision".
7. Dimension-sweep sentence reworded (always positive, largest at small sites).
8. Real bibliography added (15 entries) and A1/A2/A3 given concrete sentences.
9. Self-reference correction demoted from contribution #1 to a methodology
   subsection; contributions now lead with the diagnostics finding.

## Your task
1. Confirm each checklist item is genuinely applied (read paper/main.tex).
2. Final verdict: is this submission-plausible as a measurement paper, yes/no?
3. Residual must-fix items, if any (text only).
4. One-line summary you would give a PC member.
