"""S42 - batch-composition missingness channel on the real 1000G dosage skeleton.

s20 planted a per-marker call-rate excess for members (a hand-set Delta); s40's
M3 reduced that to a single realistic point. This gives the closed-form
sensitivity of the channel to a batch parameter: the real 1000G skeleton
supplies per-marker frequencies and population labels, and the call patterns
come from a generative model whose member-correlation is an explicit,
interpretable batch parameter.

  - panel: real 1000G dosage skeleton (three super-populations), the same
    3-way disjoint split as every other experiment.
  - batch: individuals are assigned processing batches keyed on their real 1000G
    `pop` label (1000g_panel.txt), so batch composition is population-structured.
  - quality: each individual's per-marker call probability is (1 - base_missing)
    plus a member-heavy-batch bonus, multiplied by a per-marker batch profile
    N(1, sigma_batch); member-heavy batches occur with fraction rho of the
    populations. rho therefore controls how membership-correlated the batch
    composition is: rho=0 -> no member-heavy batch (null, no membership signal);
    rho=1 -> every population's batch is member-heavy.
  - attack: the study releases its per-SNP call-rate table (fit on member batches
    + reference); a target's own called pattern is scored under the study table
    vs the reference table (the s20 M3 likelihood ratio).
  - controls: rho=0 null; a held-out self-reference control (the released table
    excludes the scored members, s20 holdout); a batch-by-super_pop control
    (coarser composition).

Sweep: rho x base-missing x sigma_batch.
Output: results/raw/s42_missingness_realistic.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import auroc
from real_multi import get_real_panel_chr, SUPERPOP

OUT = "/data/lvvyongji1/GWAS/results/raw"
PANEL = "/data/lvvyongji1/GWAS/data/external/1000g_panel.txt"
BATCH_BONUS = 0.03          # call-probability bonus for a member-heavy batch
N_SEEDS = 5                 # seeds per sweep point


def pop_labels():
    ids, pop, sup = [], [], []
    with open(PANEL) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 3:
                ids.append(p[0]); pop.append(p[1]); sup.append(p[2])
    return ids, np.array(pop), np.array(sup)


def build_split(chrom=22, per=90, seed=0):
    """Real panel + 3-way disjoint split; each selected individual's population
    label is recovered through the panel order the loader preserves."""
    G, sup, pos = get_real_panel_chr(chrom, ld_window=1500, max_snps=20000)
    ids, pops, sups = pop_labels()
    key = {s: p for s, p in zip(ids, pops)}
    order_ids = [s for s, k in zip(ids, sups) if k in SUPERPOP]
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for p in SUPERPOP:
        idx = rng.permutation(np.where(sup == p)[0])
        q = min(per, idx.size // 3)
        si.append(idx[:q]); ri.append(idx[q:2 * q]); ni.append(idx[2 * q:3 * q])
    si, ri, ni = np.concatenate(si), np.concatenate(ri), np.concatenate(ni)
    lab = lambda idx: np.array([key[order_ids[i]] for i in idx])
    return dict(G=G, sup=sup, si=si, ri=ri, ni=ni,
                pop_si=lab(si), pop_ri=lab(ri), pop_ni=lab(ni),
                sup_si=sup[si], sup_ri=sup[ri], sup_ni=sup[ni])


def batch_heavy(pop_ids, rho, seed):
    """Which individuals sit in a member-heavy batch: a fraction rho of the
    populations are flagged, by real population label."""
    pops = np.unique(pop_ids).tolist()
    rng = np.random.default_rng(seed)
    rng.shuffle(pops)
    k = int(round(rho * len(pops)))
    heavy = set(pops[:k])
    return np.array([p in heavy for p in pop_ids])


def miss_attack(cM, cN, cR, holdout=True, seed=7):
    """Released call-rate table from member batches + reference; score members vs
    non-members by the s20 M3 likelihood ratio. Also returns the candidate-only
    baseline: the AUROC of each individual's OWN total call count, which uses NO
    released table and NO reference table (it is information the candidate
    already has about itself). The incremental exposure of the release is the
    table attack minus this baseline."""
    if holdout:
        tr = np.random.default_rng(seed).random(len(cM)) < 0.5
        study_call = np.concatenate([cM[tr], cR]).mean(0)
        mem = cM[~tr]
    else:
        study_call = np.concatenate([cM, cR]).mean(0)
        mem = cM
    ref_call = cR.mean(0)
    eps = 1e-3

    def ll(called, table):
        return (called * np.log(table + eps) + (~called) * np.log(1 - table + eps)).sum(1)
    dm = ll(mem, study_call) - ll(mem, ref_call)
    dn = ll(cN, study_call) - ll(cN, ref_call)
    y = np.concatenate([np.ones(len(mem)), np.zeros(len(cN))])
    table_auc = float(auroc(np.concatenate([dm, dn]), y))
    # candidate-only: own called-marker count (higher for member-heavy batches)
    cm_only = mem.sum(1).astype(float)
    cn_only = cN.sum(1).astype(float)
    cand_auc = float(auroc(np.concatenate([cm_only, cn_only]), y))
    return table_auc, cand_auc


def patterns(nind, heavy, base_missing, m, rng, q_profile, bonus):
    """Per-individual called pattern. The per-marker quality profile q_profile is
    SHARED across panels (it is the same array's per-marker quality); only the
    per-individual batch effect (heavy) differs. A per-panel independent profile
    would create a spurious channel (a panel's own profile enters the study
    table), so sharing it is required for the null to be at chance."""
    q = np.clip((1.0 - base_missing) + bonus * heavy.astype(float), 0.5, 1.0)
    qm = np.clip(q[:, None] * q_profile[None, :], 0.3, 1.0)
    return rng.random((nind, m)) < qm


def run_cell(chrom, rho, base_missing=0.03, sigma_batch=0.02, seed=0,
             batch_by="pop", holdout=True, bonus=BATCH_BONUS):
    d = build_split(chrom=chrom, per=90, seed=seed)
    G = d["G"]
    if batch_by == "pop":
        km, kr, kn = d["pop_si"], d["pop_ri"], d["pop_ni"]
    else:
        km, kr, kn = d["sup_si"], d["sup_ri"], d["sup_ni"]
    hm = batch_heavy(km, rho, seed + 5)
    # reference and non-member panels are processed in ordinary (non-member)
    # batches: no member-heavy bonus. This is what makes the released study
    # call-rate table differ from the reference table.
    hr = np.zeros(len(kr), bool)
    hn = np.zeros(len(kn), bool)
    rng = np.random.default_rng(seed + 1000)
    m = G.shape[1]
    q_profile = np.clip(1.0 + rng.normal(0, sigma_batch, size=m), 0.3, 1.0)
    cM = patterns(len(km), hm, base_missing, m, rng, q_profile, bonus)
    cN = patterns(len(kn), hn, base_missing, m, rng, q_profile, bonus)
    cR = patterns(len(kr), hr, base_missing, m, rng, q_profile, bonus)
    return miss_attack(cM, cN, cR, holdout=holdout)


def sweep(name, values, **kw):
    out = {}
    for v in values:
        pairs = [run_cell(22, **{name: v}, seed=s, **kw) for s in range(N_SEEDS)]
        tab = [p[0] for p in pairs]; cand = [p[1] for p in pairs]
        out[f"{name}{v}"] = [float(np.mean(tab)), float(np.std(tab))]
        out[f"{name}{v}_candonly"] = [float(np.mean(cand)), float(np.std(cand))]
        print(name, v, out[f"{name}{v}"], "candonly", out[f"{name}{v}_candonly"], flush=True)
    return out


def run():
    t0 = time.time()
    res = dict(grid={}, controls={}, chrom=22)
    res["grid"]["rho"] = sweep("rho", (0.0, 0.25, 0.5, 0.75, 1.0))
    res["grid"]["base_missing"] = sweep("base_missing", (0.01, 0.02, 0.03, 0.05), rho=1.0)
    res["grid"]["sigma_batch"] = sweep("sigma_batch", (0.0, 0.01, 0.02, 0.05, 0.10), rho=1.0)
    res["grid"]["bonus"] = sweep("bonus", (0.0, 0.005, 0.01, 0.02, 0.03), rho=1.0)
    v = [run_cell(22, 1.0, seed=s, batch_by="super_pop") for s in range(N_SEEDS)]
    res["controls"]["batch_by_superpop"] = [float(np.mean([p[0] for p in v])), float(np.std([p[0] for p in v]))]
    v = [run_cell(22, 1.0, seed=s, holdout=False) for s in range(N_SEEDS)]
    res["controls"]["no_holdout"] = [float(np.mean([p[0] for p in v])), float(np.std([p[0] for p in v]))]
    print("controls", json.dumps(res["controls"]), flush=True)
    res["n_seeds"] = N_SEEDS
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s42_missingness_realistic.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
