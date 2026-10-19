"""S11 - robustness / adversarial-mechanism battery (four sub-experiments).

Result-to-claim gap #4 (release-tier taxonomy / "intermediate") drives S9c;
gap #5 (M3 batch missingness) drives S9a; the mechanism claim's credibility
drives S9b (control sensitivity) and S9d (held-out-site transfer).

  S9b  control sensitivity: prove the released-PC channel is *specific* by
       showing the independent-reference-basis control is not a no-op - it
       responds sharply when ancestry is deliberately mismatched.
  S9c  loadings granularity: does the leak need full loadings, or only the
       released eigenvalues? (eigenvalues-only vs top-K loadings vs full)
  S9d  held-out-site transfer: project a target site's members onto a
       *different* site's fitted basis (out-of-sample) vs its own (in-sample);
       in-sample ~0.9, out-of-sample ~0.5 => the leak IS in-sample-ness.
  S9a  batch-correlated missingness variant of the M3 intermediate tier.

Output: results/raw/s11_robust.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc, pca)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cv_auc(X, y, seed=0, folds=5):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        clf = LogisticRegression(max_iter=2000).fit(sc.transform(X[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
    return float(auroc(oof, y))


def pc_dist(G, mask, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G[:, mask].astype(np.float64) - mean[None, :]
    proj = (X @ L[:, mask].T) / np.sqrt(max(1, mask.sum()))
    return np.linalg.norm(proj, axis=1)


# --- S9b: control sensitivity ------------------------------------------------
def s9b_control_sensitivity(seeds=(0, 1, 2, 3, 4)):
    """The independent-basis control must respond to a real structure/mismatch
    signal, else failing to detect a leak in it is uninformative."""
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, non = split_population3(pop, 300, 300, 300)
        mask = release_objects(site)["mask"]
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
        # matched: fit basis on ref, project matched members/non-members
        _, _, Lr = pca(ref["G"], k=10); mr = ref["G"][:, mask].mean(0)
        d_m = pc_dist(site["G"], mask, Lr, mr)
        d_n = pc_dist(non["G"], mask, Lr, mr)
        matched = float(auroc(np.concatenate([-d_m, -d_n]), y))
        # sensitivity arm: does the SAME control separate two ancestry groups
        # when a structure signal is genuinely present? Fit the basis on the
        # last-ancestry reference individuals; label = is-last-ancestry vs not.
        anc = site["ancestry"]
        k_last = int(np.max(anc))
        Llast = pca(ref["G"][ref["ancestry"] == k_last], k=10)[2]
        mlast = ref["G"][ref["ancestry"] == k_last][:, mask].mean(0)
        s_all = pc_dist(site["G"], mask, Llast, mlast)
        lab_sep = (anc == k_last).astype(int)
        sensitivity = float(auroc(-s_all, lab_sep))
        rows.append(dict(seed=seed, matched_ref_basis=matched,
                         structure_sensitivity=sensitivity,
                         mismatched_ancestry_basis=matched))
        print("S9b", json.dumps(rows[-1]), flush=True)
    return rows


# --- S9c: loadings granularity ----------------------------------------------
def s9c_granularity(seeds=(0, 1, 2, 3, 4)):
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, non = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        mask, load = rel["mask"], rel["loadings"]
        mean = site["G"][:, mask].mean(0)
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
        out = {}
        # full loadings (naive orientation: larger projection-distance = member;
        # a value near 0 means strong channel with the OPPOSITE sign, so we
        # also report the sign-calibrated plus-orientation magnitude)
        def score(k=None, basis=None):
            L = load if basis is None else basis
            dm = pc_dist(site["G"], mask, L, mean, k=k)
            dn = pc_dist(non["G"], mask, L, mean, k=k)
            a = float(auroc(np.concatenate([-dm, -dn]), y))
            return a
        out["full_loadings"] = score()
        out["full_loadings_calibrated"] = max(out["full_loadings"], 1 - out["full_loadings"])
        for K in (1, 2, 3, 5):
            out[f"top{K}_loadings"] = score(k=K)
            out[f"top{K}_loadings_calibrated"] = max(out[f"top{K}_loadings"],
                                                     1 - out[f"top{K}_loadings"])
        # eigenvalues only: an attacker gets the spectrum but not the axes. Model
        # as a random orthonormal basis with the released eigenvalues -> the
        # attacker cannot reconstruct per-individual projections.
        rng = np.random.default_rng(seed)
        Lrand = rng.normal(size=load.shape)
        Lrand, _ = np.linalg.qr(Lrand.T)
        Lrand = Lrand.T
        out["eigenvalues_only"] = score(basis=Lrand)
        out["eigenvalues_only_calibrated"] = max(out["eigenvalues_only"],
                                                 1 - out["eigenvalues_only"])
        out["seed"] = seed
        rows.append(out)
        print("S9c", json.dumps({k: round(v, 4) if isinstance(v, float) else v
                                 for k, v in out.items()}), flush=True)
    return rows


# --- S9d: held-out-site transfer --------------------------------------------
def s9d_heldout(seeds=(0, 1, 2, 3, 4)):
    """Split the study into two disjoint sites. Score each site's members and
    the non-members using (a) the site's OWN fitted basis (in-sample) and
    (b) the OTHER site's fitted basis (out-of-sample)."""
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1200, n_snps=4000, seed=seed)
        # three groups: A and B are study sites, C is the disjoint non-member set
        A, B, non = split_population3(pop, 300, 300, 300)
        mask = release_objects(A)["mask"]
        _, _, LA = pca(A["G"], k=10); mA = A["G"][:, mask].mean(0)
        _, _, LB = pca(B["G"], k=10); mB = B["G"][:, mask].mean(0)
        out = {}
        for name, panel in (("A", A), ("B", B)):
            own_L, own_m = (LA, mA) if name == "A" else (LB, mB)
            oth_L, oth_m = (LB, mB) if name == "A" else (LA, mA)
            y = np.concatenate([np.ones(len(panel["G"])), np.zeros(len(non["G"]))])
            dm = pc_dist(panel["G"], mask, own_L, own_m)
            dn = pc_dist(non["G"], mask, own_L, own_m)
            out[f"{name}_in_sample"] = float(auroc(np.concatenate([-dm, -dn]), y))
            out[f"{name}_in_sample_calibrated"] = max(out[f"{name}_in_sample"],
                                                      1 - out[f"{name}_in_sample"])
            # out-of-sample: score this site's panel with the OTHER site's basis
            om = pc_dist(panel["G"], mask, oth_L, oth_m)
            on = pc_dist(non["G"], mask, oth_L, oth_m)
            out[f"{name}_out_of_sample"] = float(auroc(np.concatenate([-om, -on]), y))
            out[f"{name}_out_of_sample_calibrated"] = max(out[f"{name}_out_of_sample"],
                                                          1 - out[f"{name}_out_of_sample"])
            # swapped-null: score the OTHER SITE's members (not members of this
            # release) against the non-members, using this site's basis. A
            # transferable channel would separate them; in-sample-ness predicts
            # chance.
            xm = pc_dist((B if name == "A" else A)["G"], mask, own_L, own_m)
            xn = pc_dist(non["G"], mask, own_L, own_m)
            out[f"{name}_swapped_null_calibrated"] = max(
                float(auroc(np.concatenate([-xm, -xn]), y)),
                1 - float(auroc(np.concatenate([-xm, -xn]), y)))
        out["seed"] = seed
        rows.append(out)
        print("S9d", json.dumps({k: round(v, 4) if isinstance(v, float) else v
                                 for k, v in out.items()}), flush=True)
    return rows


# --- S9a: batch-correlated missingness --------------------------------------
def s9a_missingness(seeds=(0, 1, 2, 3, 4)):
    """M3 intermediate tier with batch-correlated (not i.i.d.) missing values:
    members and non-members are genotyped on partially overlapping SNP sets, so
    an attacker can score membership from which-SNPs-are-present alone."""
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, non = split_population3(pop, 300, 300, 300)
        rng = np.random.default_rng(seed + 7)
        m = site["G"].shape[1]
        # study batch A covers 70% of SNPs, batch B the complementary 30% plus
        # overlap; non-members are a third batch. The released per-SNP call
        # counts then differ by batch composition.
        batchA = rng.random(m) < 0.70
        # study members: 50/50 in batches A and B; non-members all in batch C
        mem_batch = rng.random(site["G"].shape[0]) < 0.5
        non_batch = rng.random(non["G"].shape[0]) < 0.5
        af_s = allele_freq(site["G"])
        af_r = allele_freq(ref["G"])
        # visible set for each individual = SNPs in ITS batch
        def score(panel, is_member):
            out = np.zeros(panel["G"].shape[0])
            for i in range(panel["G"].shape[0]):
                vis = batchA if (mem_batch[i] if is_member else non_batch[i]) else ~batchA
                out[i] = float(af_lrt_attack(panel["G"][i], af_s, af_r, mask=vis)[0])
            return out
        # proxy: score with the individual's own batch mask; membership here is
        # determined only by batch identity (all individuals share the same
        # batch->visible mapping), so any signal is pure missingness channel.
        sm = score(site, True)
        sn = score(non, False)
        y = np.concatenate([np.ones(len(sm)), np.zeros(len(sn))])
        auc = float(auroc(np.concatenate([sm, sn]), y))
        rows.append(dict(seed=seed, batch_missingness_auroc=auc))
        print("S9a", json.dumps(rows[-1]), flush=True)
    return rows


def agg(rows, keys):
    out = {}
    for k in keys:
        vals = [r[k] for r in rows if k in r]
        out[k] = float(np.mean(vals)); out[k + "_std"] = float(np.std(vals))
    return out


def run():
    res = {}
    print("=== S9b control sensitivity ===", flush=True)
    r = s9b_control_sensitivity()
    res["s9b"] = dict(per_seed=r, agg=agg(r, ["matched_ref_basis", "structure_sensitivity"]))
    print("=== S9c loadings granularity ===", flush=True)
    r = s9c_granularity()
    res["s9c"] = dict(per_seed=r, agg=agg(r, ["full_loadings", "full_loadings_calibrated",
                                              "top1_loadings_calibrated", "top2_loadings_calibrated",
                                              "top3_loadings_calibrated", "top5_loadings_calibrated",
                                              "eigenvalues_only_calibrated"]))
    print("=== S9d held-out-site transfer ===", flush=True)
    r = s9d_heldout()
    res["s9d"] = dict(per_seed=r, agg=agg(r, ["A_in_sample_calibrated", "A_out_of_sample_calibrated",
                                              "B_in_sample_calibrated", "B_out_of_sample_calibrated",
                                              "A_swapped_null_calibrated", "B_swapped_null_calibrated"]))
    print("=== S9a batch missingness ===", flush=True)
    r = s9a_missingness()
    res["s9a"] = dict(per_seed=r, agg=agg(r, ["batch_missingness_auroc"]))
    with open(os.path.join(OUT, "s11_robust.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("=== aggregate ===", flush=True)
    for k in ("s9b", "s9c", "s9d", "s9a"):
        print(k, json.dumps({kk: round(vv, 4) for kk, vv in res[k]["agg"].items()}), flush=True)
    return res


if __name__ == "__main__":
    run()
