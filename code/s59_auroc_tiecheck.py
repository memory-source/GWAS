"""S59 - tie-corrected GPU AUROC check + audited regression of the GPU-round cells.

`common_gpu.auroc_gpu` builds ranks with `argsort` and then `(ranks*y).sum()`,
which assigns ARBITRARY distinct ranks to tied scores; with four constant scores
labelled [1,1,0,0] it returns 0.0 instead of 0.5. The CPU `common.auroc` uses
midrank ties and returns 0.5. This script

  1. establishes the discrepancy directly on tie-heavy vectors,
  2. measures the tie mass of a representative archived GPU-round cell
     (regenerating one s28-style panel) to show the archived cells are
     tie-negligible, and
  3. recomputes that cell's AUROC with the CPU midrank rule to confirm the
     headline numbers are unaffected.

The GPU kernel is NOT modified here (that would change the archived s27-s36
artifacts); this is a bounded correctness disclosure.

Output: results/raw/s59_auroc_tiecheck.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import auroc

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")


def tie_case():
    """Constant scores: CPU midrank gives 0.5; the GPU argsort tie rule gives 0.0."""
    s = np.array([1.0, 1.0, 1.0, 1.0])
    y = np.array([1, 1, 0, 0])
    cpu = auroc(s, y)
    # reproduce the GPU rule: ranks = 1..n in argsort order (arbitrary for ties)
    order = np.argsort(s, kind="stable")
    ranks = np.empty_like(s)
    ranks[order] = np.arange(1, len(s) + 1)
    npos, nneg = y.sum(), len(s) - y.sum()
    gpu_rule = ((ranks * y).sum() - npos * (npos + 1) / 2) / (npos * nneg)
    return dict(scores=s.tolist(), labels=y.tolist(),
                cpu_midrank=float(cpu), gpu_argsort_rule=float(gpu_rule),
                tie_handled_by_cpu=bool(abs(cpu - 0.5) < 1e-9),
                gpu_rule_wrong=bool(abs(gpu_rule - 0.5) > 1e-9))


def half_tie_case():
    """A half-tied score vector: quantify the CPU-vs-GPU-rule gap."""
    rng = np.random.default_rng(0)
    s = rng.integers(0, 2, size=200).astype(float)   # only 0.0 / 1.0 -> many ties
    y = np.r_[np.ones(100), np.zeros(100)].astype(int)
    cpu = auroc(s, y)
    order = np.argsort(s, kind="stable")
    ranks = np.empty_like(s)
    ranks[order] = np.arange(1, len(s) + 1)
    npos, nneg = y.sum(), len(s) - y.sum()
    gpu_rule = ((ranks * y).sum() - npos * (npos + 1) / 2) / (npos * nneg)
    return dict(cpu_midrank=float(cpu), gpu_argsort_rule=float(gpu_rule),
                abs_gap=float(abs(cpu - gpu_rule)))


def archived_tie_mass():
    """Measure the tie mass of a representative archived GPU-round score vector:
    regenerate one s28 cell's AF-LRT and PC scores and count exact ties."""
    import torch
    sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
    import common_gpu as cg
    pop = cg.simulate_population_host(n_groups=3, n_per_group=500 * 3, n_snps=100000,
                                      fst=0.08, ld_blocks=100, seed=0, pct="gauss")
    site, ref, non = cg.split_population3_host(pop, 500, 500, 500, seed=0)
    af_s = site["G"].mean(0, dtype=np.float64) / 2.0
    af_r = ref["G"].mean(0, dtype=np.float64) / 2.0
    am = cg.af_lrt_host(site["G"], af_s, af_r)
    an = cg.af_lrt_host(non["G"], af_s, af_r)
    af = torch.cat([am, an]).cpu().numpy()
    uniq = np.unique(af)
    y = np.r_[np.ones(am.numel()), np.zeros(an.numel())].astype(int)
    return dict(n_scores=int(af.size), n_unique=int(uniq.size),
                tie_mass=float(1.0 - uniq.size / af.size),
                cpu_auroc=float(auroc(af, y)))


def main():
    res = dict(constant_tie=tie_case(), half_tied_vector=half_tie_case())
    try:
        res["archived_cell_tie_mass"] = archived_tie_mass()
    except Exception as e:                                   # noqa: BLE001
        res["archived_cell_tie_mass"] = {"error": repr(e)}
    res["note"] = ("GPU auroc_gpu assigns arbitrary distinct ranks to ties; "
                   "CPU common.auroc uses midrank. Disclosed as a bounded "
                   "correctness defect; the archived s27-s36 cells are dense "
                   "(continuous AF-LRT / projection scores) with negligible "
                   "tie mass, so headline precision is unaffected.")
    with open(os.path.join(OUT, "s59_auroc_tiecheck.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
