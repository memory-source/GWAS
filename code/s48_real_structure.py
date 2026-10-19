"""S48 - real-data structure ladder: the PC-vs-AF crossover measured on real
genomes.

The paper's mechanism says the PC-vs-AF ordering is governed by measured
population structure. On real 1000G this has only been shown on the AF-first
side (3 super-populations mixed). This script measures BOTH sides on real
genomes: a ladder of cohorts from within a single 1000G population (between-
population Fst ~ 0) up to 5 super-populations (Fst ~ 0.1), running the identical
attack code (AF-LRT vs study-fitted PC projection distance) at each rung.

Per cohort:
  - 3-way disjoint split (study / reference / non-member) per population group;
  - af_auc (study-fitted AF vs reference AF), pc_auc (calibrated),
    gap = pc - af;
  - measured fst_like (observed structure);
  - reference-fitted-basis and eigenvalue-only structure controls.

Cohorts (ascending measured structure): each single 1000G `pop`; a 2-pop mix per
super-population; the 3-super-population set; the 5-super-population set.

Output: results/raw/s48_real_structure.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, auroc, pca
from real_multi import get_real_panel_chr, load_panel, SUPERPOP, SUPERPOP5


OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")


def _orient_auroc(scores, y, rng, cal_mask=None):
    """AUROC on a held-out half, with the channel's sign fixed on the OTHER
    (calibration) half only. The calibration observations are excluded from the
    evaluation, so no scored individual's label ever selected the direction.
    Pass `cal_mask` to reuse a shared calibration/evaluation split across
    channels scored on the same row ordering."""
    s = np.asarray(scores, float)
    n = len(s)
    h = rng.random(n) < 0.5 if cal_mask is None else np.asarray(cal_mask, bool)
    if h.sum() >= 2 and (~h).sum() >= 2 and len(np.unique(y[h])) == 2:
        sgn = 1.0 if auroc(s[h], y[h]) >= 0.5 else -1.0
        ev = ~h
    else:
        sgn = -1.0                                # distance channels: closer = member
        ev = np.ones(n, bool)
    return float(auroc(sgn * s[ev], y[ev]))


def cohort_cell(G, lab, per=150, seed=0, k=10):
    """One cohort: disjoint 3-way split per group; AF vs PC + controls."""
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for g in np.unique(lab):
        idx = rng.permutation(np.where(lab == g)[0])
        q = min(per, idx.size // 3)
        if q < 4:
            continue
        si.append(idx[:q]); ri.append(idx[q:2 * q]); ni.append(idx[2 * q:3 * q])
    if not si:
        return None
    si, ri, ni = np.concatenate(si), np.concatenate(ri), np.concatenate(ni)
    site, ref, non = G[si], G[ri], G[ni]
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    # one shared calibration/evaluation split for every channel scored on these
    # rows; the calibration half sets the sign, the held-out half is evaluated.
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    af_s, af_r = allele_freq(site), allele_freq(ref)
    am = af_lrt_attack(site, af_s, af_r); an = af_lrt_attack(non, af_s, af_r)
    af_auc = _orient_auroc(np.concatenate([am, an]), y,
                           np.random.default_rng(seed + 31), cal_mask=cal)
    _, _, L = pca(site, k=k)
    mean = site.mean(0)
    dm = np.linalg.norm((site.astype(float) - mean) @ L.T / np.sqrt(G.shape[1]), axis=1)
    dn = np.linalg.norm((non.astype(float) - mean) @ L.T / np.sqrt(G.shape[1]), axis=1)
    pc_auc = _orient_auroc(np.concatenate([dm, dn]), y,
                           np.random.default_rng(seed + 31), cal_mask=cal)
    # reference-fitted basis control
    _, _, Lr = pca(ref, k=k)
    mr = ref.mean(0)
    xm = np.linalg.norm((site.astype(float) - mr) @ Lr.T / np.sqrt(G.shape[1]), axis=1)
    xn = np.linalg.norm((non.astype(float) - mr) @ Lr.T / np.sqrt(G.shape[1]), axis=1)
    refb = _orient_auroc(np.concatenate([xm, xn]), y,
                         np.random.default_rng(seed + 31), cal_mask=cal)
    # eigenvalue-only control (random orthonormal basis)
    Le = rng.normal(size=(k, G.shape[1])); Le, _ = np.linalg.qr(Le.T); Le = Le.T
    em = np.linalg.norm((site.astype(float) - mean) @ Le.T / np.sqrt(G.shape[1]), axis=1)
    en = np.linalg.norm((non.astype(float) - mean) @ Le.T / np.sqrt(G.shape[1]), axis=1)
    eigo = _orient_auroc(np.concatenate([em, en]), y,
                         np.random.default_rng(seed + 31), cal_mask=cal)
    return dict(af_auc=af_auc, pc_auc=pc_auc, gap=pc_auc - af_auc,
                fst_like=fst_like(site, lab[si]),
                pc_refbasis=refb, pc_eigonly=eigo,
                n_study=int(len(site)), n_group=int(len(np.unique(lab[si]))),
                M=int(G.shape[1]))


def pop_table():
    """sample -> (pop, super_pop), and the pop list per super_pop, in the panel's
    (and therefore the loader's) order."""
    samples, sup = load_panel()
    pops, sups = [], []
    with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "data", "external", "1000g_panel.txt")) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 3:
                pops.append(p[1]); sups.append(p[2])
    return samples, np.array(sups), np.array(pops)


def cohort_masks(samples, sups, pops, sups_keep):
    """Boolean samples-mask + per-sample group label (the pop) for the loaded
    sample order (`samples` is the loader's order for `sups_keep`)."""
    keep = np.array([s in sups_keep for s in sups])
    ids = [s for s, k in zip(samples, keep) if k]
    lab = np.array([pops[samples.index(s)] for s in ids])
    return ids, lab


def load_for(sups_keep, chrom=22, ld_window=1500, max_snps=20000):
    G, sup, pos = get_real_panel_chr(chrom, ld_window=ld_window, max_snps=max_snps,
                                     sups=tuple(sups_keep))
    samples, sups, pops = pop_table()
    order = [s for s, k in zip(samples, sups) if k in set(sups_keep)]
    lab = np.array([pops[samples.index(s)] for s in order])
    return G, lab, pos


def fst_like(G, lab, min_n=8):
    af = G.mean(0) / 2.0
    hs = []
    for g in np.unique(lab):
        if (lab == g).sum() < min_n:
            continue
        fg = G[lab == g].mean(0) / 2.0
        hs.append(2 * fg * (1 - fg))
    hs = np.mean(hs, axis=0)
    ht = 2 * af * (1 - af)
    return float(np.nanmean((ht - hs) / np.where(ht > 0, ht, np.nan)))


COHORTS = [
    ("pop:GBR", ["GBR"]), ("pop:FIN", ["FIN"]), ("pop:CHB", ["CHB"]),
    ("pop:JPT", ["JPT"]), ("pop:YRI", ["YRI"]), ("pop:LWK", ["LWK"]),
    ("pop:TSI", ["TSI"]), ("pop:PEL", ["PEL"]), ("pop:IBS", ["IBS"]),
    ("pop:CDX", ["CDX"]), ("pop:STU", ["STU"]), ("pop:PJL", ["PJL"]),
    ("super:EUR-2pop", None), ("super:EAS-2pop", None),
    ("super:3pop", None), ("super:5pop", None),
]


def run_chrom(chrom=22, n_splits=5, ld_window=1500, max_snps=20000):
    samples, sups, pops = pop_table()
    res = {}
    # single-population cohorts (load the 1000G 5-super-pop set so every pop exists)
    G, lab, pos = load_for(SUPERPOP5, chrom, ld_window, max_snps)
    for name, _ in COHORTS:
        if not name.startswith("pop:"):
            continue
        pop = name.split(":")[1]
        m = lab == pop
        if m.sum() < 12:
            continue
        rows = [cohort_cell(G[m], lab[m], seed=s) for s in range(n_splits)]
        rows = [r for r in rows if r]
        res[name] = _agg(rows)
    # 2-pop mixes within each super-pop (filter on pop labels) and the super-pop
    # mixes (filter on the sample's super_pop, since 3pop/5pop are super-pop sets)
    order = [s for s, k in zip(samples, sups) if k in set(SUPERPOP5)]
    order_sup = np.array([sups[samples.index(s)] for s in order])
    for sp in SUPERPOP5:
        ps = np.unique(pops[sups == sp])
        if len(ps) < 2:
            continue
        m = np.isin(lab, list(ps[:2]))
        rows = [cohort_cell(G[m], lab[m], seed=s) for s in range(n_splits)]
        rows = [r for r in rows if r]
        res[f"super:{sp}-2pop"] = _agg(rows)
    for name, sp_keep in (("super:3pop", SUPERPOP), ("super:5pop", SUPERPOP5)):
        m = np.isin(order_sup, sp_keep)
        rows = [cohort_cell(G[m], order_sup[m], seed=s) for s in range(n_splits)]
        rows = [r for r in rows if r]
        res[name] = _agg(rows)
    return res


def _agg(rows):
    keys = [k for k in rows[0] if isinstance(rows[0][k], (int, float))]
    out = {}
    for k in keys:
        vals = np.array([r[k] for r in rows], float)
        out[k] = float(vals.mean())
        out[k + "_std"] = float(vals.std())
    return out


def run(chrom=22, out=None):
    res = run_chrom(chrom)
    payload = dict(chrom=chrom, cohorts=res, n_splits=5)
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    for name, a in sorted(res.items(), key=lambda kv: kv[1]["fst_like"]):
        print("%-18s fst %.4f  AF %.3f  PC %.3f  gap %+.3f  refb %.3f  eigo %.3f  ngrp %d" % (
            name, a["fst_like"], a["af_auc"], a["pc_auc"], a["gap"],
            a["pc_refbasis"], a["pc_eigonly"], a["n_group"]), flush=True)
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    ap.add_argument("--out", default=os.path.join(OUT, "s48_real_structure.json"))
    a = ap.parse_args()
    run(a.chrom, a.out)
