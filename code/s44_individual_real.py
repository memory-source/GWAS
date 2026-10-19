"""S44 - per-individual privacy metric on real and semi-synthetic data.

Round-11's per-individual metric (s38) was computed only on the parametric
simulator. This re-expresses the same estimators on two real-structure panels:

  - real 1000G chr22 (real_multi.get_real_panel_chr), the same 3-way disjoint
    split as s16/three_way_split;
  - a semi-synthetic cohort (semi_synthetic.generate_cohort at nsp=2000), which
    carries real 1000G haplotype structure but is large enough for the 3-way
    split.

Channels: AF-LRT at full precision and 1dp, the per-super-population (sitewise)
AF, the study-fitted PC projection, and the eigenvalue-only / reference-fitted
controls. The metric is s38's `metrics_from_scores`, now reporting both the
`resolution` FNR floor (1/n_pos) and the confidence-corrected `rule3` floor
(3/n_pos) so the f-DP eps_lb is quotable where the empirical TPR reaches 1.

Output: results/raw/s44_individual_real.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, pca, auroc
from s38_privacy_metric import metrics_from_scores, pc_dist


OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")


def _agg(acc):
    out = {}
    for name, rows in acc.items():
        out[name] = {}
        for k in rows[0]:
            vals = [r[k] for r in rows if isinstance(r[k], (int, float, bool))]
            out[name][k] = float(np.mean(vals))
        for k in rows[0]:
            if k.endswith("_valid"):
                out[name][k] = bool(all(r[k] for r in rows))
    return out


def metrics_cell(G_site, G_ref, G_non, seed, sup_labels=None):
    """The s38 channel battery on one already-split panel triple."""
    mask = np.ones(G_site.shape[1], bool)
    af_s_full = allele_freq(G_site)
    af_r = allele_freq(G_ref)
    y = np.concatenate([np.ones(len(G_site)), np.zeros(len(G_non))])
    out = {}
    for name, af in (("af_full", af_s_full),
                     ("af_1dp", np.round(af_s_full, 1))):
        sm = af_lrt_attack(G_site, af, af_r)
        sn = af_lrt_attack(G_non, af, af_r)
        out[name] = metrics_from_scores(np.concatenate([sm, sn]), y)
    # per-super-population (sitewise) AF
    if sup_labels is not None:
        sm = np.zeros(len(G_site)); sn = np.zeros(len(G_non))
        for g in np.unique(sup_labels):
            afk = allele_freq(G_site[sup_labels == g])
            rf = allele_freq(G_ref[sup_labels == g])
            sm[sup_labels == g] = af_lrt_attack(G_site[sup_labels == g], afk, rf)
            sn[sup_labels == g] = af_lrt_attack(G_non[sup_labels == g], afk, rf)
        out["af_sitewise"] = metrics_from_scores(np.concatenate([sm, sn]), y)
    # study-fitted PC
    _, _, L = pca(G_site, k=10)
    mean = G_site.mean(0)
    dm = pc_dist(G_site, L, mean); dn = pc_dist(G_non, L, mean)
    rng = np.random.default_rng(seed + 7717)
    tsm = rng.random(len(dm)) < 0.5; tsn = rng.random(len(dn)) < 0.5
    ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
    sgn = 1.0 if auroc(-np.concatenate([dm[tsm], dn[tsn]]), ytr) >= 0.5 else -1.0
    out["pc_loadings"] = metrics_from_scores(sgn * (-np.concatenate([dm, dn])), y)
    # eigenvalue-only control (random orthonormal basis)
    Lr = rng.normal(size=L.shape); Lr, _ = np.linalg.qr(Lr.T); Lr = Lr.T
    dmr = pc_dist(G_site, Lr, mean); dnr = pc_dist(G_non, Lr, mean)
    out["pc_eigonly_control"] = metrics_from_scores(sgn * (-np.concatenate([dmr, dnr])), y)
    # reference-fitted basis control
    _, _, Lref = pca(G_ref, k=10)
    mr = G_ref.mean(0)
    dmx = pc_dist(G_site, Lref, mr); dnx = pc_dist(G_non, Lref, mr)
    out["pc_refbasis_control"] = metrics_from_scores(-np.concatenate([dmx, dnx]), y)
    return out


def real_chr(chrom=22, n_splits=8, per=150):
    from real_multi import get_real_panel_chr, three_way_split
    G, sup, pos = get_real_panel_chr(chrom)
    acc = {}
    for s in range(n_splits):
        si, ri, ni = three_way_split(sup, per=per, seed=s)
        cells = metrics_cell(G[si], G[ri], G[ni], s, sup_labels=sup[si])
        for name, m in cells.items():
            acc.setdefault(name, []).append(m)
        print("chr%d split%d done" % (chrom, s), flush=True)
    return dict(chrom=chrom, n_snps=int(G.shape[1]), agg=_agg(acc), n_splits=n_splits)


def cohort_chr(chrom=22, nsp=2000, rho=1e-7, M=8000, seeds=(0, 1, 2)):
    """One generated cohort per seed; the per-seed channel battery is pooled."""
    acc = {}
    for seed in seeds:
        c = cohort_one(chrom, nsp, rho, seed, M)
        for name, m in c.items():
            acc.setdefault(name, []).append(m)
        print("cohort seed%d done" % seed, flush=True)
    return dict(agg=_agg(acc), nsp=nsp, M=M, seeds=list(seeds))


def cohort_one(chrom=22, nsp=2000, rho=1e-7, seed=0, M=8000):
    from semi_synthetic import generate_cohort
    c = generate_cohort(chrom, n_per_group=nsp, rho=rho, seed=seed, stride_to=M)
    G, anc = c["G"], c["ancestry"]
    # 3-way disjoint per group
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for k in range(3):
        idx = rng.permutation(np.where(anc == k)[0])
        q = min(nsp, idx.size // 3)
        si.append(idx[:q]); ri.append(idx[q:2 * q]); ni.append(idx[2 * q:3 * q])
    si, ri, ni = np.concatenate(si), np.concatenate(ri), np.concatenate(ni)
    return metrics_cell(G[si], G[ri], G[ni], seed, sup_labels=anc[si])


def run():
    import time
    path = os.path.join(OUT, "s44_individual_real.json")
    res = {}
    if os.path.exists(path):
        try:
            res = json.load(open(path))
        except Exception:
            res = {}
    res["space"] = os.environ.get("S44_SPACE", "both")
    t0 = time.time()
    if res["space"] in ("both", "real"):
        res["real_chr22"] = real_chr(22)
        print("real chr22", json.dumps({k: round(v["auroc"], 3)
              for k, v in res["real_chr22"]["agg"].items()}), flush=True)
    if res["space"] in ("both", "cohort"):
        res["cohort_nsp2000"] = cohort_chr(22, nsp=2000)
        print("cohort nsp2000", json.dumps({k: round(v["auroc"], 3)
              for k, v in res["cohort_nsp2000"]["agg"].items() if isinstance(v, dict)}), flush=True)
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    print("done", res["elapsed_s"], "s", flush=True)
    return res


if __name__ == "__main__":
    run()
