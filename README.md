# GWAS — released diagnostics, membership leakage, and reveal governance

Measurement study of what a multi-site genome-wide association study (GWAS) releases.
The project has three parts:

1. **Passive leakage measurement.** Given a fixed release (allele frequencies, the
   study-fitted PC loading matrix, eigenvalues, missingness, haplotype tables), measure
   how much *membership* information each released object carries, using a controlled
   protocol in which members, the attacker's reference panel, and the scored non-members
   are disjoint.
2. **Selection-input spectral tomography.** A cooperating analysis partner changes its
   own input across joint-PCA re-fits while keeping every per-SNP genotype count fixed,
   reads only the truncated final output, and accumulates membership information beyond
   the strongest passive adversary built from the same releases.
3. **Governance layer.** A pre-registered reveal manifest, an enforced reveal egress
   monitor, and a hash-chained cross-run ledger, built and stress-tested against careless
   authoring and a native-code boundary.

Every number in `paper/` traces to one JSON file in `results/raw/`.

## Headline results

Passive (paper in `paper/main.tex`, claims in `CLAIMS_FROM_RESULTS.md`):

- The released study-fitted PC **loading matrix** is a stronger membership channel than
  the unconditioned allele-frequency summary in the small-(n × M) regime: projection
  distance 0.88 AUROC vs 0.84. The mechanism is the in-sample-ness of the fitted
  spectrum (eigenvalue-only at chance; refit-on-fresh-panel at chance).
- A conditioned-AF attack (subtract known calibration members from the released summary)
  reaches 0.941, above the full-study PC channel, so the AF ordering is stated
  explicitly instead of claiming universal PC dominance.
- Scoring the attacker's own reference individuals as non-members inflates AUROC by
  about 0.13; the paper's controlled protocol removes that.

Active (round 22, claims in `CLAIMS_FROM_RESULTS_ARB_20261017.md`):

- The tomography identity `Y_t = V_t Λ_t − K_t V_t = C_H V_t` holds on real data
  (residual 3.8e-11 in float64, 6.2e-6 through a real PLINK 2.0 export read-back).
- An AF-preserving block swap leaves per-SNP histograms, allele frequency and variance
  exactly unchanged; only the individual correspondence moves.
- Stacking q = 8 outputs gives a membership increment over the strongest of four passive
  constructions (conditioned AF, conditioned single release, their label-free fusion and
  their calibration-weighted fusion), evaluated as a held-out oracle:
  chr22 AUC gap +0.0066 (EUR), +0.0656 (3 populations), +0.0037 (EAS);
  chr1 +0.0130 / +0.0443 / +0.0048; positive interval in all cohorts.
- In a **matched printed-output** comparison (`s72`) — the same q = 8 re-fits through
  PLINK 2.0 `--pca 10 allele-wts approx`, read back and compared against the passive
  bar built from the same printed outputs — the AUC gap is positive in every cell
  (+0.0054 EUR chr22, +0.0082 3pop chr22, +0.0100 EUR chr1, +0.0060 3pop chr1), but the
  conditional split interval excludes zero only for chr1 EUR. The matched TPR@5%FPR
  reading is not established.

- Bounded precision: the reproducible-direction count under simulated quantization falls
  from 78 to 30 (6 decimal places) and 11 (4 decimal places) on chr22 EUR; that count is
  a heuristic, not a usable-rank theorem.

## What is not established

- Independent-subject replication (chr1 reuses the same 1000G subjects).
- Biobank-scale behavior of PLINK's randomized `approx` solver (at n = 240 it agrees with
  the exact solver, so the solver-robustness result is a lower bound).
- Printed-output cross-population superiority: positive mean, interval spans zero.
- A clean null Type-I-error rate for the scientific-cost arm: the reported quantity is
  the off-causal discovery fraction, and LD proxies of causal markers can contribute.

## Repository layout

| path | contents |
|---|---|
| `code/` | 91 Python modules: `s0`–`s72` experiment scripts plus helpers (`common.py`, `real_panel.py`, `governance.py`, `enforce.py`, `common_gpu.py`, `make_figures.py`) |
| `results/raw/` | one JSON per experiment cell; every quoted number comes from here |
| `paper/` | `main.tex` + `main.pdf` + `figures/` |
| `refine-logs/` | experiment plans, results write-ups, tracker |
| `review-stage/` | cross-model review logs, reviewer memory, round state |
| `research-plan/` | research directions and the deep-search specification |
| `idea-stage/` | idea discovery report and the research contract |
| `EXPERIMENT_AUDIT*.md/json` | integrity audits (one per round) |
| `proposal_experiment_plan.md` | the original Chinese proposal |
| `data/` | data provenance (`data/README.md`); VCFs are not tracked |

## Reproducing

Environment:

```
mamba create -p ./env python=3.10 numpy=1.26 scipy scikit-learn pandas matplotlib seaborn
./env/bin/pip install msprime==1.3.1 tskit==0.5.6 scikit-allel
```

Data: download the 1000G Phase-3 VCFs and the HapMap genetic maps as described in
`data/README.md`, then place them under `data/external/`.

PLINK 2.0: `code/s63_plink_export.py` and `code/s72_printed_matched.py` expect a
`plink2` binary at `scratch/plink2bin/plink2`. Any PLINK 2.0 build works; the printed
allele weights and eigenvalues are read from `.eigenvec.allele` / `.eigenval`.

Run a script:

```
OMP_NUM_THREADS=4 ./env/bin/python code/s67_joint_passive.py --chrom 22
```

Panels are cached to `results/cache/` after the first parse (the 1000G VCFs take a few
minutes per chromosome). All scripts write to `results/raw/`.

## Review and audit trail

- `EXPERIMENT_AUDIT_ARB.md` / `EXPERIMENT_AUDIT_ARB.json` — round-22 integrity audit.
- `review-stage/ARB_REVIEW_20261017.md` — the four-round auto-review log that ended at
  "ready" for the narrowed claim.
- `review-stage/DEEP_REVIEW_20261016.md` — the ten-round review of the deep-search
  direction.
- `review-stage/REVIEW_STATE*.json` — machine-readable round state.
- `CLAIMS_FROM_RESULTS*.md` — claim maps with per-item verdicts, including what is
  explicitly *not* supported.

## Citation

If you use this code or measurement protocol, cite the paper draft in `paper/main.tex`.
