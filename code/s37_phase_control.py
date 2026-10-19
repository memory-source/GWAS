"""Relatedness control with the ACTUAL copied haplotype phase.

The s37/s44 cohort controls reported a dosage parsimony-phase proxy for
relatedness, which understates the true copied-haplotype sharing. This probe
regenerates the same cohorts with generate_cohort(..., return_hap=True) and
reports the longest-shared-run relatedness on the ACTUAL mosaic haplotypes
(cohort vs its own donors), alongside the dosage proxy, so the matched-
relatedness claim can be stated on the real phase.

Output: results/raw/s37_phase_control.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
import semi_synthetic as ss
from real_phased import SUPERPOP

OUT = "/data/lvvyongji1/GWAS/results/raw"
CELLS = [(150, 8000), (2000, 8000)]


def one(chrom, nsp, M, seed=0):
    d = ss.load_donors(chrom, stride_to=M)
    c = ss.generate_cohort(chrom, n_per_group=nsp, rho=1e-7, seed=seed,
                           stride_to=M, return_hap=True)
    H = np.vstack([c["H0"], c["H1"]])
    don_lib = np.concatenate([d["lib"][p] for p in SUPERPOP], axis=0)
    Dd = np.concatenate([d["dosage"][p] for p in SUPERPOP], axis=0)
    rel_actual_c = ss.relatedness(H, seed=seed)
    rel_actual_d = ss.relatedness(don_lib, seed=seed)
    rel_proxy_c = ss.relatedness(ss.parsimony_phase(c["G"]), seed=seed)
    rel_proxy_d = ss.relatedness(ss.parsimony_phase(Dd), seed=seed)
    return dict(nsp=nsp, M=M,
                actual_cohort=rel_actual_c, actual_donor=rel_actual_d,
                proxy_cohort=rel_proxy_c, proxy_donor=rel_proxy_d,
                actual_ratio=float(rel_actual_c[0] / rel_actual_d[0]),
                proxy_ratio=float(rel_proxy_c[0] / rel_proxy_d[0]))


def run():
    out = {}
    for (nsp, M) in CELLS:
        out[f"nsp{nsp}_m{M}"] = one(22, nsp, M)
        print(f"nsp{nsp}_m{M}", json.dumps({k: (round(v, 3) if isinstance(v, float)
              else [round(x, 2) for x in v] if isinstance(v, list) else v)
              for k, v in out[f"nsp{nsp}_m{M}"].items()}), flush=True)
    with open(os.path.join(OUT, "s37_phase_control.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    run()
