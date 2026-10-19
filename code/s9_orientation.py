"""S9 - PC-channel orientation study, real 1000G chr22 + matched-operating-point
simulation.

Motivation (result-to-claim gap #1): the simulation's study-fitted
PC-distance channel is strongly *directional* - in-sample study members sit
farther from the released centroid than disjoint non-members, which is the
+0.10 increment. On the real panel the same naive distance is anti-predictive
(0.347). This experiment asks, at a MATCHED operating point (450/panel,
~1.4k SNPs), whether:

  (a) which orientation the PC channel takes, and
  (b) whether the sim's injected LD blocks (ld_blocks>0) are what make the
      channel so strongly directional,

and reports the three orientations an attacker could hold:
  naive    = -projection-norm  (members-farther hypothesis)
  flipped  = +projection-norm
  oracle   = max(naive, flipped)  (attack-side sign calibration, uses labels)

Output: results/raw/s9_orientation.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, allele_freq,
                    af_lrt_attack, auroc, pca)
from real_panel import get_real_panel, three_way_split

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cv_auc(X, y, seed=0, folds=5):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        clf = LogisticRegression(max_iter=2000)
        clf.fit(sc.transform(X[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
    return float(auroc(oof, y))


def proj_norm(G, mask, loadings, mean):
    X = G[:, mask].astype(np.float64) - mean[None, :]
    proj = (X @ loadings[:, mask].T) / np.sqrt(max(1, mask.sum()))
    return np.linalg.norm(proj, axis=1)


def orientation_block(site_G, non_G, mask, loadings, mean, af_s, af_r, y, seed):
    """Return the per-channel AUROCs at one operating point."""
    d_m = proj_norm(site_G, mask, loadings, mean)
    d_n = proj_norm(non_G, mask, loadings, mean)
    # naive = "members sit farther from the centroid" (score = +||proj||);
    # flipped = the opposite orientation (score = -||proj||).
    naive = float(auroc(np.concatenate([d_m, d_n]), y))
    flipped = float(auroc(np.concatenate([-d_m, -d_n]), y))
    oracle = max(naive, flipped)
    af_m = af_lrt_attack(site_G, af_s, af_r, mask=mask)
    af_n = af_lrt_attack(non_G, af_s, af_r, mask=mask)
    af_only = float(auroc(np.concatenate([af_m, af_n]), y))
    # combined: AF + *calibrated* PC (orientation fixed by the labeled split,
    # as an attacker with a train split would do)
    orientation = "farther" if naive >= flipped else "closer"
    signed_m = d_m if orientation == "farther" else -d_m
    signed_n = d_n if orientation == "farther" else -d_n
    comb = cv_auc(np.column_stack([np.concatenate([af_m, af_n]),
                                   np.concatenate([signed_m, signed_n])]), y, seed)
    return dict(af_only=af_only, pc_naive=naive, pc_flipped=flipped,
                pc_oracle=oracle, af_plus_calibrated_pc=comb,
                pc_orientation=orientation)


def real_arm(split_seeds=tuple(range(20))):
    G, sup, pos = get_real_panel()          # (n_samples, ~1363)
    mask = np.ones(G.shape[1], bool)
    rows = []
    for s in split_seeds:
        sidx, ridx, nidx = three_way_split(sup, per=150, seed=s)
        site, ref, non = G[sidx], G[ridx], G[nidx]
        _, _, load_s = pca(site, k=10)
        mean_s = site.mean(0)
        af_s = allele_freq(site)
        af_r = allele_freq(ref)
        y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
        r = orientation_block(site, non, mask, load_s, mean_s, af_s, af_r, y, s)
        r["seed"] = s
        rows.append(r)
        print("real seed", s, json.dumps(r), flush=True)
    return rows


def sim_arm(ld_blocks_list=(0, 50, 100, 200), n_seeds=20, n_groups=3,
            n_per_group=450, n_snps=1363):
    out = {}
    for lb in ld_blocks_list:
        rows = []
        for seed in range(n_seeds):
            pop = simulate_population(n_groups=n_groups, n_per_group=n_per_group,
                                      n_snps=n_snps, ld_blocks=lb, seed=seed)
            site, ref, non = split_population3(pop, 150, 150, 150)
            rel_mask = np.ones(site["G"].shape[1], bool)
            _, _, load_s = pca(site["G"], k=10)
            mean_s = site["G"].mean(0)
            af_s = allele_freq(site["G"])
            af_r = allele_freq(ref["G"])
            y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
            r = orientation_block(site["G"], non["G"], rel_mask, load_s, mean_s,
                                  af_s, af_r, y, seed)
            r["seed"] = seed
            rows.append(r)
        agg = {}
        for k in ("af_only", "pc_naive", "pc_flipped", "pc_oracle",
                  "af_plus_calibrated_pc"):
            agg[k] = float(np.mean([x[k] for x in rows]))
            agg[k + "_std"] = float(np.std([x[k] for x in rows]))
        agg["orientation_farther_frac"] = float(
            np.mean([x["pc_orientation"] == "farther" for x in rows]))
        out[f"ld_blocks={lb}"] = dict(agg=agg, per_seed=rows)
        print("sim ld_blocks", lb, json.dumps({k: round(v, 4) for k, v in agg.items()}),
              flush=True)
    return out


def run():
    os.makedirs(OUT, exist_ok=True)
    print("=== real arm ===", flush=True)
    real_rows = real_arm()
    real_agg = {}
    for k in ("af_only", "pc_naive", "pc_flipped", "pc_oracle",
              "af_plus_calibrated_pc"):
        real_agg[k] = float(np.mean([x[k] for x in real_rows]))
        real_agg[k + "_std"] = float(np.std([x[k] for x in real_rows]))
    real_agg["orientation_farther_frac"] = float(
        np.mean([x["pc_orientation"] == "farther" for x in real_rows]))
    print("real agg", json.dumps({k: round(v, 4) for k, v in real_agg.items()}), flush=True)

    print("=== sim arm (matched op point) ===", flush=True)
    sim = sim_arm()

    res = dict(real=dict(agg=real_agg, per_seed=real_rows),
               sim_matched=sim,
               meta=dict(n_per_panel=450, real_n_snps=int(get_real_panel()[0].shape[1]),
                         real_split_seeds=20, sim_seeds=20, sim_n_snps=1363))
    with open(os.path.join(OUT, "s9_orientation.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
