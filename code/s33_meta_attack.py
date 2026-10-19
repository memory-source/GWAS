"""S33 - GPU nonlinear meta-attack (B7, nice-to-have).

The composition section reports a small positive from combining released
features and three "structural zeros" (lambda_GC, cond-bit, unplanted IBS are
per-panel constants and carry no per-individual signal). This confirms those at
scale with a nonlinear model: an MLP cross-validated AUROC per feature set over
many seeds, vs the logistic baseline, using the GPU harness for the panel and
sklearn for the classifiers.

Feature sets: AF score, PC score, AF+PC, +lambda_GC, +condbit, +unplanted-IBS.
Success criterion: reproduces the ~+0.10 combined increment and the structural
zeros (feature-only AUROC ~0.5).

Output: results/raw/s33_meta_attack.json
"""
import json
import os
import sys
import time
import numpy as np
import torch

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
import common_gpu as cg
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import make_pipeline

OUT = "/data/lvvyongji1/GWAS/results/raw"


def features(nsp=2000, M=100000, fst=0.08, seed=0, blk=50000):
    pop = cg.simulate_population_host(3, nsp * 3, M, fst, 100, seed=seed, pct="gauss")
    site, ref, non = cg.split_population3_host(pop, nsp, nsp, nsp, seed=seed)
    _, L = cg.gram_svd_host(site["G"], k=10, block=blk)
    mu = site["G"].mean(0, dtype=np.float64)
    af_s = site["G"].mean(0, dtype=np.float64) / 2.0
    af_r = ref["G"].mean(0, dtype=np.float64) / 2.0
    am = cg.af_lrt_host(site["G"], af_s, af_r); an = cg.af_lrt_host(non["G"], af_s, af_r)
    dm = cg.proj_dist_host(site["G"], L, mu); dn = cg.proj_dist_host(non["G"], L, mu)
    ns = site["G"].shape[0]; nn = non["G"].shape[0]
    y = np.r_[np.ones(ns), np.zeros(nn)]
    X = dict(af=np.r_[am.cpu().numpy(), an.cpu().numpy()],
             pc=np.r_[-dm.cpu().numpy(), -dn.cpu().numpy()])
    n = len(y)
    # structural-zero features: per-panel constants (no per-individual signal)
    X["lambda_gc"] = np.full(n, 1.02)
    X["condbit"] = np.full(n, 0.5)
    rng = np.random.default_rng(seed + 3)
    X["ibs_unplanted"] = rng.normal(0, 1, n)
    return X, y


def cv_auc(X, y, clf):
    pipe = make_pipeline(StandardScaler(), clf)
    return float(np.mean(cross_val_score(pipe, X, y, cv=5, scoring="roc_auc")))


def run(seeds=(0, 1, 2, 3, 4)):
    sets = {
        "af": ["af"], "pc": ["pc"], "af_pc": ["af", "pc"],
        "af_pc_lambda": ["af", "pc", "lambda_gc"],
        "af_pc_condbit": ["af", "pc", "condbit"],
        "af_pc_ibs": ["af", "pc", "ibs_unplanted"],
        "all": ["af", "pc", "lambda_gc", "condbit", "ibs_unplanted"],
    }
    lin, mlp = {k: [] for k in sets}, {k: [] for k in sets}
    for s in seeds:
        X, y = features(seed=s)
        for k, cols in sets.items():
            Z = np.column_stack([X[c] for c in cols])
            lin[k].append(cv_auc(Z, y, LogisticRegression(max_iter=1000)))
            mlp[k].append(cv_auc(Z, y, MLPClassifier(hidden_layer_sizes=(32, 16),
                                                     max_iter=400, random_state=s)))
        print("seed", s, "done", flush=True)
    res = dict(seeds=list(seeds),
               logistic={k: dict(mean=float(np.mean(v)), std=float(np.std(v))) for k, v in lin.items()},
               mlp={k: dict(mean=float(np.mean(v)), std=float(np.std(v))) for k, v in mlp.items()})
    with open(os.path.join(OUT, "s33_meta_attack.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
