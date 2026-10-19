"""S1c - frontier (disjoint negatives, 5 seeds, per-ancestor), efficient:
simulate each seed's panel once and reuse it across all tiers.

Output: results/raw/s1c_frontier.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects, allele_freq,
    af_lrt_attack, anc_only_attack, auroc, tpr_at_fpr, pca,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")
F_TIERS = ["F2_topk", "F4_round1", "F3_round2", "F5_full", "F5_sitewise"]


def released_af(rel, tier):
    if tier == "F4_round1":
        return np.round(rel["af"], 1)
    if tier == "F3_round2":
        return np.round(rel["af"], 2)
    return rel["af"]


def visible(rel, tier, topk=200):
    base = rel["mask"].copy()
    if tier == "F2_topk":
        keep = np.zeros_like(base)
        tested = np.where(base)[0]
        keep[tested[np.argsort(rel["pvals"][tested])][:topk]] = True
        return keep
    return base


def bootstrap_ci(scores, labels, fpr=0.01, n_boot=300, seed=0):
    rng = np.random.default_rng(seed)
    pos = np.where(labels == 1)[0]; neg = np.where(labels == 0)[0]
    vals = []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(pos, pos.size, True), rng.choice(neg, neg.size, True)])
        vals.append(tpr_at_fpr(scores[idx], labels[idx], fpr))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def run(seeds=(0, 1, 2, 3, 4)):
    # simulate once per seed
    panels = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, nonmem = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        panels.append((site, ref, nonmem, rel))
        print("simulated seed", seed, flush=True)

    res = {}
    for tier in F_TIERS:
        aurocs, tpr1, lo_l, hi_l, anc, study_pc = [], [], [], [], [], []
        per_anc = {0: [], 1: [], 2: []}
        for (site, ref, nonmem, rel) in panels:
            vis = visible(rel, tier)
            af_s = released_af(rel, tier)
            af_r = allele_freq(ref["G"])
            if tier == "F5_sitewise":
                sm = np.zeros(len(site["G"])); sn = np.zeros(len(nonmem["G"]))
                for k in rel["af_by_group"]:
                    afk = rel["af_by_group"][k]
                    sm[site["ancestry"] == k] = af_lrt_attack(
                        site["G"][site["ancestry"] == k], afk,
                        allele_freq(ref["G"][ref["ancestry"] == k]), mask=rel["mask"])
                    sn[nonmem["ancestry"] == k] = af_lrt_attack(
                        nonmem["G"][nonmem["ancestry"] == k], afk,
                        allele_freq(ref["G"][ref["ancestry"] == k]), mask=rel["mask"])
            else:
                sm = af_lrt_attack(site["G"], af_s, af_r, mask=vis)
                sn = af_lrt_attack(nonmem["G"], af_s, af_r, mask=vis)
            lab = np.concatenate([np.ones(len(sm)), np.zeros(len(sn))])
            sc = np.concatenate([sm, sn])
            aurocs.append(auroc(sc, lab))
            tpr1.append(tpr_at_fpr(sc, lab, 0.01))
            lo, hi = bootstrap_ci(sc, lab, seed=0)
            lo_l.append(lo); hi_l.append(hi)
            for k in range(3):
                sk = np.concatenate([sm[site["ancestry"] == k], sn[nonmem["ancestry"] == k]])
                lk = np.concatenate([np.ones((site["ancestry"] == k).sum()),
                                     np.zeros((nonmem["ancestry"] == k).sum())])
                per_anc[k].append(auroc(sk, lk))
            # Ancestry/structure confound control: fit the PC basis on an
            # INDEPENDENT panel (the attacker reference), project both target
            # panels, and score proximity to the reference centroid. Because
            # members and non-members are ancestry-matched and out-of-sample
            # w.r.t. this basis, this isolates the pure structure channel with
            # no study-fitting leak and no AF information.
            _, _, load_ref = pca(ref["G"], k=10)
            mean_ref = ref["G"][:, rel["mask"]].mean(0)
            av_m = anc_only_attack(site["G"], load_ref, mean_ref, rel["mask"])
            av_n = anc_only_attack(nonmem["G"], load_ref, mean_ref, rel["mask"])
            anc.append(auroc(np.concatenate([av_m, av_n]), lab))
            # Study-fitted PC-distance attack (the S2c leak, for comparison)
            _, _, load_s = pca(site["G"], k=10)
            mean_s = site["G"][:, rel["mask"]].mean(0)
            sv_m = anc_only_attack(site["G"], load_s, mean_s, rel["mask"])
            sv_n = anc_only_attack(nonmem["G"], load_s, mean_s, rel["mask"])
            study_pc.append(auroc(np.concatenate([sv_m, sv_n]), lab))
        res[tier] = dict(auroc=float(np.mean(aurocs)), auroc_std=float(np.std(aurocs)),
                         tpr_at_fpr1=float(np.mean(tpr1)),
                         tpr_ci=[float(np.mean(lo_l)), float(np.mean(hi_l))],
                         ancestry_only_auroc=float(np.mean(anc)),
                         study_pc_distance_auroc=float(np.mean(study_pc)),
                         per_ancestor={int(k): float(np.mean(v)) for k, v in per_anc.items()})
        print(tier, json.dumps(res[tier]), flush=True)
    with open(os.path.join(OUT, "s1c_frontier.json"), "w") as f:
        json.dump(dict(f_tiers=F_TIERS, results=res, seeds=list(seeds)), f, indent=2)
    return res


if __name__ == "__main__":
    run()
