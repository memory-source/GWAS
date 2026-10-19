"""S46 - second bundle manifest + uncoarsened single channel.

s40 measured one realistic manifest (pooled AF@2dp, per-site AF@2dp, top-10
eig, lambda_GC@2dp, call-rate table, top-50 pvals) against the strongest single
object (per-site AF@2dp). Two robustness gaps remain:

  1. the result could be specific to that one manifest's object set / precisions;
  2. the single-object comparator was coarsened, so an "oracle max" against an
     UNCOARSENED per-site AF might erase the excess.

This script runs the same joint out-of-fold logistic attacker under a SECOND
manifest (per-site AF at FULL precision + top-10 eigenvalues + lambda_GC@2dp +
the call-rate table; pooled AF and p-values omitted), and computes the strongest
single object both coarsened (as in s40) and uncoarsened. Leave-one-out
attribution and the structural-zero noise control are kept.

Output: results/raw/s46_bundle_manifest2.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)
from governance import apply_precision
from s40_bundle import make_call_table, m3_from_table

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")

MANIFEST2 = {
    "af_by_group": {"kind": "full"},
    "eig": {"kind": "topk", "k": 10},
    "lambda_gc": {"kind": "decimals", "d": 2},
    "missing": {"kind": "decimals", "d": 2},
}


def build_bundle2(rel):
    b = {}
    for q, spec in MANIFEST2.items():
        if spec["kind"] == "full":
            b[q] = rel[q] if q != "af_by_group" else {k: v for k, v in rel[q].items()}
        else:
            b[q] = apply_precision(rel[q], spec)
    b["mask"] = rel["mask"]
    return b


def features2(site, ref, non, rel, m3, coarsen_site):
    """Feature matrix under manifest 2. `coarsen_site` rounds the per-site AF to
    2dp (the s40 comparator) when True, else uses the full-precision per-site AF."""
    b = build_bundle2(rel)
    mask = rel["mask"]
    f_as_m = np.full(len(site["G"]), -1e9)
    f_as_n = np.full(len(non["G"]), -1e9)
    for k in b["af_by_group"]:
        afk = b["af_by_group"][k]
        if coarsen_site:
            afk = np.round(afk, 2)
        smk = af_lrt_attack(site["G"][site["ancestry"] == k], afk,
                            allele_freq(ref["G"][ref["ancestry"] == k]), mask=mask)
        snk = af_lrt_attack(non["G"][non["ancestry"] == k], afk,
                            allele_freq(ref["G"][ref["ancestry"] == k]), mask=mask)
        f_as_m[site["ancestry"] == k] = smk
        f_as_n[non["ancestry"] == k] = snk
    dm, dn = m3
    return (dict(af_site=f_as_m, m3=dm),
            dict(af_site=f_as_n, m3=dn))


def joint_auc(fm, fn, seed, cols):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    Xm = np.column_stack([fm[c] for c in cols])
    Xn = np.column_stack([fn[c] for c in cols])
    X = np.vstack([Xm, Xn]); y = np.concatenate([np.ones(len(Xm)), np.zeros(len(Xn))])
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for a, b in skf.split(X, y):
        sc = StandardScaler().fit(X[a]); clf = LogisticRegression(max_iter=2000)
        clf.fit(sc.transform(X[a]), y[a]); oof[b] = clf.predict_proba(sc.transform(X[b]))[:, 1]
    return float(auroc(oof, y))


COLS = ["af_site", "m3"]


def run(seeds=(0, 1, 2, 3, 4)):
    rows = []
    for s in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=s)
        site, ref, non = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        # released call-rate table (the missingness object); attacker consumes
        # its coarsened @2dp emit, matching MANIFEST2.
        call_table, cm, cn, cr, tr = make_call_table(site, ref, non, s)
        rel["missing"] = call_table
        m3 = m3_from_table(np.round(call_table, 2), cm, cn, cr)
        # coarsened single-object comparator (s40 style) and uncoarsened
        fm_c, fn_c = features2(site, ref, non, rel, m3, coarsen_site=True)
        fm_u, fn_u = features2(site, ref, non, rel, m3, coarsen_site=False)
        # candidate-only comparator column (own call count; no released table)
        fm_c["m3_cand"] = cm.sum(1).astype(float); fn_c["m3_cand"] = cn.sum(1).astype(float)
        fm_u["m3_cand"] = cm.sum(1).astype(float); fn_u["m3_cand"] = cn.sum(1).astype(float)
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
        def orient(a):
            return max(float(a), 1 - float(a))
        single_coarse = orient(auroc(np.concatenate([fm_c["af_site"], fn_c["af_site"]]), y))
        single_uncoarse = orient(auroc(np.concatenate([fm_u["af_site"], fn_u["af_site"]]), y))
        bundle2 = joint_auc(fm_u, fn_u, s, COLS)
        loo = {c: joint_auc(fm_u, fn_u, s, [x for x in COLS if x != c]) for c in COLS}
        # structural-zero control: random per-individual column added
        rng = np.random.default_rng(s + 1234)
        fm2 = dict(fm_u, noise=rng.normal(size=len(site["G"])))
        fn2 = dict(fn_u, noise=rng.normal(size=len(non["G"])))
        with_noise = joint_auc(fm2, fn2, s, ["af_site", "noise"])
        single_plus_noise = joint_auc(fm2, fn2, s, ["af_site"])
        cand_bundle = joint_auc(fm_u, fn_u, s, ["af_site", "m3_cand"])
        cand_plus_table = joint_auc(fm_u, fn_u, s, ["af_site", "m3_cand", "m3"])
        rows.append(dict(seed=s, single_coarse=single_coarse, single_uncoarse=single_uncoarse,
                         bundle2=bundle2, excess_vs_coarse=bundle2 - single_coarse,
                         excess_vs_uncoarse=bundle2 - single_uncoarse,
                         loo=loo, single_plus_noise=single_plus_noise,
                         af_site_plus_noise=with_noise,
                         cand_bundle=cand_bundle, cand_plus_table=cand_plus_table,
                         m3_increment=cand_plus_table - cand_bundle))
        print("seed", s, "bundle2 %.3f single_coarse %.3f single_uncoarse %.3f excess_u %+.3f" % (
            bundle2, single_coarse, single_uncoarse, bundle2 - single_uncoarse), flush=True)
    agg = dict(
        bundle2=float(np.mean([r["bundle2"] for r in rows])),
        single_coarse=float(np.mean([r["single_coarse"] for r in rows])),
        single_uncoarse=float(np.mean([r["single_uncoarse"] for r in rows])),
        excess_vs_coarse=float(np.mean([r["excess_vs_coarse"] for r in rows])),
        excess_vs_uncoarse=float(np.mean([r["excess_vs_uncoarse"] for r in rows])),
        excess_uncoarse_frac_gt0=float(np.mean([r["excess_vs_uncoarse"] > 0 for r in rows])),
        loo={c: float(np.mean([r["loo"][c] for r in rows])) for c in COLS},
        single_plus_noise=float(np.mean([r["single_plus_noise"] for r in rows])),
        af_site_plus_noise=float(np.mean([r["af_site_plus_noise"] for r in rows])),
        cand_bundle=float(np.mean([r["cand_bundle"] for r in rows])),
        cand_plus_table=float(np.mean([r["cand_plus_table"] for r in rows])),
        m3_increment=float(np.mean([r["m3_increment"] for r in rows])),
        m3_increment_frac_gt0=float(np.mean([r["m3_increment"] > 0 for r in rows])))
    res = dict(per_seed=rows, agg=agg, n_seeds=len(seeds), manifest=MANIFEST2, cols=COLS)
    with open(os.path.join(OUT, "s46_bundle_manifest2.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps({k: round(v, 3) for k, v in agg.items() if not isinstance(v, dict)}), flush=True)
    return res


if __name__ == "__main__":
    run()
