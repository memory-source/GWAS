"""S40 - joint release-bundle leakage.

s23 measures ONE governed object in isolation; s24 a declared SEQUENCE with a
cumulative logistic. Neither answers whether a realistic multi-object bundle
leaks more than its parts. This script measures the FULL aggregate of a
realistic bundle with a joint attacker that consumes every released object, and
reports the excess over the strongest single object.

Bundle (from the existing release_objects output, coarsened by apply_precision):
  pooled AF@2dp, per-site AF@2dp, top-10 eigenvalues, lambda_GC@2dp,
  per-SNP call-rate table@2dp (the s20 M3 channel), top-50 p-values. The
  call-rate table the study releases IS the object the attacker consumes.

Controls:
  - per-panel-constant guard: eig + lambda_GC + pvals carry no per-individual
    signal, so they cannot enter a per-individual classifier; the guard is a
    random per-individual noise column stacked onto pooled AF, which must not
    raise the AUROC;
  - single-object governed negative: a manifest releasing only pooled AF@2dp
    must collapse the aggregate to that object's value.

Output: results/raw/s40_bundle.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)
from governance import apply_precision

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")

MANIFEST = {
    "af": {"kind": "decimals", "d": 2},
    "af_by_group": {"kind": "decimals", "d": 2},
    "eig": {"kind": "topk", "k": 10},
    "lambda_gc": {"kind": "decimals", "d": 2},
    "missing": {"kind": "decimals", "d": 2},
    "pvals": {"kind": "topk", "k": 50, "order": "small"},
}


def build_bundle(rel):
    b = {}
    for q, spec in MANIFEST.items():
        if q == "af_by_group":
            b[q] = {k: apply_precision(v, spec) for k, v in rel[q].items()}
        else:
            b[q] = apply_precision(rel[q], spec)
    b["mask"] = rel["mask"]
    return b


def make_call_table(site, ref, non, seed, mem_rate=0.031, non_rate=0.030):
    """The released per-SNP call-rate object. Members carry a small higher
    per-marker call-failure rate (a membership-correlated genotyping-batch
    effect; realistic chip missingness magnitudes); the table is fit on a TRAIN
    half of the study plus the reference. Returns the table (the released
    object), the per-individual called patterns the attacker observes, and the
    train mask `tr` so a caller can separate members who CONTRIBUTED to the
    table from members who did not (the released table is a study release, so
    every member's pattern is observable to the attacker)."""
    rng = np.random.default_rng(seed + 909)
    m = site["G"].shape[1]
    nm, nn = len(site["G"]), len(non["G"])
    tr = rng.random(nm) < 0.5
    cm = rng.random((nm, m)) > mem_rate
    cn = rng.random((nn, m)) > non_rate
    cr = rng.random((len(ref["G"]), m)) > non_rate
    table = np.concatenate([cm[tr], cr]).mean(0)
    return table, cm, cn, cr, tr


def m3_from_table(table, cm, cn, cr):
    """Attacker's membership score = log-likelihood of an individual's called
    pattern under the RELEASED study table minus under the reference table. The
    table passed here is the coarsened object the manifest actually emits; the
    attacker never sees a higher-precision surrogate."""
    eps = 1e-3
    ref_table = cr.mean(0)
    table = np.asarray(table)

    def ll(called, t):
        return (called * np.log(t + eps) + (~called) * np.log(1 - t + eps)).sum(1)

    return ll(cm, table) - ll(cm, ref_table), ll(cn, table) - ll(cn, ref_table)


def features(site, ref, non, rel, seed, callpat):
    """Per-individual feature matrix from ALL released objects. The missingness
    feature is computed from `b["missing"]` -- the coarsened table the manifest
    emits -- not from a separate higher-precision table."""
    b = build_bundle(rel)
    mask = rel["mask"]
    af_r = allele_freq(ref["G"])
    f_af_m = af_lrt_attack(site["G"], b["af"], af_r, mask=mask)
    f_af_n = af_lrt_attack(non["G"], b["af"], af_r, mask=mask)
    f_as_m = np.full(len(site["G"]), -1e9)
    f_as_n = np.full(len(non["G"]), -1e9)
    for k in b["af_by_group"]:
        afk = b["af_by_group"][k]
        smk = af_lrt_attack(site["G"][site["ancestry"] == k], afk,
                            allele_freq(ref["G"][ref["ancestry"] == k]), mask=mask)
        snk = af_lrt_attack(non["G"][non["ancestry"] == k], afk,
                            allele_freq(ref["G"][ref["ancestry"] == k]), mask=mask)
        f_as_m[site["ancestry"] == k] = smk
        f_as_n[non["ancestry"] == k] = snk
    cm, cn, cr = callpat[:3]
    dm, dn = m3_from_table(b["missing"], cm, cn, cr)
    # candidate-only comparator: each individual's own called-marker count, which
    # uses NO released table and NO reference table. The m3 feature's incremental
    # value is what it adds beyond this candidate-side information.
    cand_m = cm.sum(1).astype(float); cand_n = cn.sum(1).astype(float)
    feats_m = dict(af=f_af_m, af_site=f_as_m, m3=dm, m3_cand=cand_m)
    feats_n = dict(af=f_af_n, af_site=f_as_n, m3=dn, m3_cand=cand_n)
    return feats_m, feats_n


def joint_auc(feats_m, feats_n, seed, cols):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    Xm = np.column_stack([feats_m[c] for c in cols])
    Xn = np.column_stack([feats_n[c] for c in cols])
    X = np.vstack([Xm, Xn]); y = np.concatenate([np.ones(len(Xm)), np.zeros(len(Xn))])
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for a, b in skf.split(X, y):
        sc = StandardScaler().fit(X[a]); clf = LogisticRegression(max_iter=2000)
        clf.fit(sc.transform(X[a]), y[a]); oof[b] = clf.predict_proba(sc.transform(X[b]))[:, 1]
    return float(auroc(oof, y))


COLS_ALL = ["af", "af_site", "m3"]
COLS_STRUCT = []   # eig/lambda/pvals are per-panel constants -> no per-individual column


def run(seeds=(0, 1, 2, 3, 4)):
    rows = []
    for s in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=s)
        site, ref, non = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        # the missingness object the study releases IS the call-rate table; the
        # attacker consumes its coarsened emit, not a higher-precision surrogate.
        call_table, cm, cn, cr, tr = make_call_table(site, ref, non, s)
        rel["missing"] = call_table
        fm, fn = features(site, ref, non, rel, s, (cm, cn, cr))
        # contributor split: the released study table includes the train half of
        # members, so the attacker can score all members; report M3 AUROC among
        # the non-contributor members separately (a stricter, held-out reading).
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
        single = {c: float(auroc(np.concatenate([fm[c], fn[c]]), y)) for c in COLS_ALL}
        for c in single:
            single[c] = max(single[c], 1 - single[c])
        fm_nc = dict(fm, m3=np.where(tr, -1e9, fm["m3"]))
        m3_noncontrib = float(auroc(np.concatenate([fm_nc["m3"][~tr], fn["m3"]]),
                                    np.concatenate([np.ones((~tr).sum()), np.zeros(len(non["G"]))])))
        bundle = joint_auc(fm, fn, s, COLS_ALL)
        loo = {c: joint_auc(fm, fn, s, [x for x in COLS_ALL if x != c]) for c in COLS_ALL}
        af_only_bundle = joint_auc(fm, fn, s, ["af"])
        # candidate-information comparator: AF channels + the candidate-only call
        # count (no released table), then the same bundle plus the released-table
        # M3 score. The increment is the conditional exposure of releasing the
        # call-rate table beyond what the candidate already observes.
        cand_cols = ["af", "af_site", "m3_cand"]
        cand_bundle = joint_auc(fm, fn, s, cand_cols)
        cand_plus_table = joint_auc(fm, fn, s, cand_cols + ["m3"])
        m3_increment = cand_plus_table - cand_bundle
        # per-panel-constant guard: a random per-individual column stacked onto
        # pooled AF must not raise the AUROC (the paper's eig/lambda/pvals are
        # per-panel constants with no per-individual variation and cannot enter
        # a per-individual classifier at all).
        rng = np.random.default_rng(s + 1234)
        fm2 = dict(fm, noise=rng.normal(size=len(site["G"])))
        fn2 = dict(fn, noise=rng.normal(size=len(non["G"])))
        af_plus_noise = joint_auc(fm2, fn2, s, ["af", "noise"])
        rows.append(dict(seed=s, single=single, bundle=bundle,
                         max_single=max(single.values()),
                         argmax=max(single, key=single.get),
                         excess=bundle - max(single.values()),
                         loo=loo, af_only_bundle=af_only_bundle,
                         af_plus_noise=af_plus_noise,
                         m3_noncontrib=m3_noncontrib,
                         cand_bundle=cand_bundle,
                         cand_plus_table=cand_plus_table,
                         m3_increment=m3_increment))
        print("seed", s, "bundle %.3f max_single %.3f (%s) excess %+.3f  af+noise %.3f" % (
            bundle, max(single.values()), max(single, key=single.get),
            bundle - max(single.values()), af_plus_noise), flush=True)
    agg = dict(
        bundle=float(np.mean([r["bundle"] for r in rows])),
        max_single=float(np.mean([r["max_single"] for r in rows])),
        excess=float(np.mean([r["excess"] for r in rows])),
        excess_frac_gt0=float(np.mean([r["excess"] > 0 for r in rows])),
        single={c: float(np.mean([r["single"][c] for r in rows])) for c in COLS_ALL},
        loo={c: float(np.mean([r["loo"][c] for r in rows])) for c in COLS_ALL},
        af_only_bundle=float(np.mean([r["af_only_bundle"] for r in rows])),
        af_plus_noise=float(np.mean([r["af_plus_noise"] for r in rows])),
        m3_noncontrib=float(np.mean([r["m3_noncontrib"] for r in rows])),
        cand_bundle=float(np.mean([r["cand_bundle"] for r in rows])),
        cand_plus_table=float(np.mean([r["cand_plus_table"] for r in rows])),
        m3_increment=float(np.mean([r["m3_increment"] for r in rows])),
        m3_increment_frac_gt0=float(np.mean([r["m3_increment"] > 0 for r in rows])))
    res = dict(per_seed=rows, agg=agg, n_seeds=len(seeds),
               manifest=MANIFEST, cols=COLS_ALL)
    with open(os.path.join(OUT, "s40_bundle.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps({k: (round(v, 3) if isinstance(v, float) else v)
                             for k, v in agg.items() if not isinstance(v, dict)}), flush=True)
    return res


if __name__ == "__main__":
    run()
