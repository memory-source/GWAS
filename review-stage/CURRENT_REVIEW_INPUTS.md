# Round C10 review inputs (current state)

Changed since the C9 (round 23, thread 01a12007) review:

- code/s10_paired.py — replaced the selfref/LOO framing. Now: `af_conditioned`
  (balanced known-member residual, one frequency, disjoint negatives; 0.941),
  `pc_conditioned` (PCA REFIT on the SAME remaining panel, single predefined
  orientation; 0.964), `af_loo_oracle` (member-only leave-own-out, labeled an
  oracle ablation; 0.429). pc_conditioned beats af_conditioned on 25/25 seeds
  (paired +0.023, 95% CI [+0.019,+0.028]). Results: results/raw/s10_paired.json.
- paper/main.tex — "Paired significance and an auxiliary-information control"
  paragraph rewritten around the symmetric budget (conditioned AF vs refit PC)
  and the panel-mismatch explanation (full-study PC vs conditioned AF); abstract
  and contributions updated; eigenvalue sentence reworded; discussion
  "constructed bundle"; bundle attacker split (four-of-five folds vs half-study
  contributor set); UDP probe (fixed 127.0.0.1:9, stopped at socket creation);
  Table 1 caption "per-seed bootstrap bounds for TPR@1%"; fig5(c) y-tick shows
  study individuals.
- code/make_figures.py — fig6(d) axis "study individuals"; fig5(c) y-ticks
  show study individuals; fig8(a) per-(size,Fst) curves.
- code/s27_gpu_scale.py — `combined` joint AF+PC; shards 0/1/2/3 RE-RUN with the
  corrected estimator (results/raw/s27_gpu_scale_shard{0,1,2,3}.json); the true
  AF-only for the saturated cells in results/raw/s27_afhigh_shard5.json.

Artifacts: results/raw/s10_paired.json, s27_gpu_scale_shard{0..3}.json,
s27_afhigh_shard5.json, s39_enforcement.json, s38_privacy_metric.json.
Paper: 20 pages, 0 errors, 0 undefined refs.
