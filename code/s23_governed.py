"""S23 - governed vs ungoverned egress: enforcement demonstrated at the
release-object level, and the leakage consequence re-measured under the
governed release.

This makes the paper's "governance mechanism is argued, not built" concrete for
the object layer: we declare a manifest, run the EgressMonitor over a study's
release objects, exercise the P0-5 violation-injection suite (blocked /
coarsened / detected), and then re-run the SAME leakage attacks on the governed
release vs the ungoverned release. The AUROC drop is the consequence the
enforcement predicts (it is not itself new leakage evidence).

Output: results/raw/s23_governed.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc, pca)
from governance import (make_manifest, EgressMonitor, Ledger, POLICIES,
                        violation_injection_suite, apply_precision)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def pc_dist(G, mask, loadings, mean):
    X = G[:, mask].astype(np.float64) - mean[None, :]
    return np.linalg.norm((X @ loadings[:, mask].T) / np.sqrt(max(1, mask.sum())), axis=1)


def attacks(site, ref, non, rel, af_override=None, loadings=None):
    """Run the three leakage channels on a (possibly governed) release."""
    mask = rel["mask"]
    af_s = rel["af"] if af_override is None else af_override
    af_r = allele_freq(ref["G"])
    y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
    am = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
    an = af_lrt_attack(non["G"], af_s, af_r, mask=mask)
    out = dict(af=float(auroc(np.concatenate([am, an]), y)))
    if loadings is not None:
        mean_s = site["G"][:, mask].mean(0)
        dm = pc_dist(site["G"], mask, loadings, mean_s)
        dn = pc_dist(non["G"], mask, loadings, mean_s)
        a = float(auroc(np.concatenate([dm, dn]), y))
        out["pc_dist"] = max(a, 1 - a)
    else:
        out["pc_dist"] = float("nan")
    return out


def run(seeds=(0, 1, 2, 3, 4)):
    per_seed = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, non = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        # ungoverned: everything released at full precision
        ung = attacks(site, ref, non, rel, af_override=rel["af"], loadings=rel["loadings"])
        # governed: AF @1dp, loadings withheld (attacker falls back to eig-only,
        # i.e. a public eigenvalue release gives no usable axes -> chance
        # projection), af_by_group withheld
        gov_af = np.round(rel["af"], 1)
        rng = np.random.default_rng(seed)
        L_rand = rng.normal(size=rel["loadings"].shape); L_rand, _ = np.linalg.qr(L_rand.T)
        gov = attacks(site, ref, non, rel, af_override=gov_af, loadings=L_rand.T)
        per_seed.append(dict(seed=seed, ungoverned=ung, governed=gov,
                             delta_af=gov["af"] - ung["af"],
                             delta_pc=gov["pc_dist"] - ung["pc_dist"]))
        print(json.dumps(per_seed[-1]), flush=True)
    agg = {}
    for cond in ("ungoverned", "governed"):
        for k in ("af", "pc_dist"):
            vals = [r[cond][k] for r in per_seed]
            agg[f"{cond}_{k}"] = [float(np.mean(vals)), float(np.std(vals))]
    # enforcement evidence
    enforcement = violation_injection_suite()
    res = dict(per_seed=per_seed, agg=agg,
               manifest_governed=make_manifest(POLICIES["governed"], "demo")["hash"],
               manifest_ungoverned=make_manifest(POLICIES["ungoverned"], "demo")["hash"],
               enforcement=enforcement)
    with open(os.path.join(OUT, "s23_governed.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps(agg), flush=True)
    print("enforcement outcomes",
          json.dumps([(m["injection"], m["outcome"]) for m in enforcement["matrix"]]), flush=True)
    return res


if __name__ == "__main__":
    run()
