"""S62 - honest-sample-size sensitivity for the selection-input spectral tomography.

The deep-review flagged that the first version grew BOTH honest and partner and let
split_roles silently truncate the shadow role (540 requested > 503 EUR rows). This
version holds partner / nonmember / reference / shadow FIXED and grows ONLY the
honest cohort, using a nested design: the n=120 honest rows are a subset of the
n=180 honest rows, and every other role is identical across the two arms.

EUR chr22 has 503 rows at the 1.5 kb thinning; the base assignment uses 420 and the
extra 60 honest rows come from the remaining pool (503 - 420 = 83 >= 60).

Output: results/raw/s62_sample_size_chr22.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s61_spectral_tomography as s61

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
BASE = {"honest": 120, "partner": 120, "nonmember": 80, "reference": 60, "shadow": 40}
N_HONEST_BIG = 180


def nested_roles(groups, seed):
    """Nested split: base roles from BASE, honest grown to N_HONEST_BIG from the
    unused remainder. Returns (roles_120, roles_180) with identical partner /
    nonmember / reference / shadow and honest_120 ⊂ honest_180."""
    rng = np.random.default_rng(seed)
    roles = {r: [] for r in BASE}
    rest = {r: [] for r in BASE}
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        start = 0
        for r, cap in BASE.items():
            if start + cap > idx.size:
                raise ValueError(f"group {g} has {idx.size} rows, need {start + cap}")
            roles[r].append(idx[start:start + cap])
            start += cap
        rest.setdefault("_unused", []).append(idx[start:])
    roles = {r: np.concatenate(v) for r, v in roles.items()}
    unused = np.concatenate(rest["_unused"])
    extra = N_HONEST_BIG - BASE["honest"]
    if unused.size < extra:
        raise ValueError(f"need {extra} extra honest rows, only {unused.size} unused")
    rng2 = np.random.default_rng(seed + 500)
    add = rng2.permutation(unused)[:extra]
    big = dict(roles)
    big["honest"] = np.concatenate([roles["honest"], add])
    return roles, big


def _boot_ci(diff):
    rng = np.random.default_rng(7)
    boot = diff[rng.integers(0, len(diff), size=(10000, len(diff)))].mean(1)
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]


def run(chrom=22, seeds=5):
    G, lab, _ = s48.load_for(("EUR",), chrom, 1500, 20000)
    groups = np.array(["EUR"] * len(G))
    if groups.size != 503:
        raise ValueError(f"expected 503 EUR chr{chrom} rows, got {groups.size}")
    res = {"chrom": chrom, "seeds": seeds, "marker_count": int(G.shape[1]), "sizes": {}}
    for tag, which in (("per120", 0), ("per180", 1)):
        rows = []
        for s in range(seeds):
            roles = nested_roles(groups, s)[which]
            rows.append(s61.cell(G, groups, s, roles=roles))
            if s == 0:
                first_roles = roles
        afc = np.array([r["af_conditioned_auc"] for r in rows])
        spc = np.array([r["spec_conditioned_auc"] for r in rows])
        strong = np.maximum(afc, spc)
        block = {"n_honest": int(len(first_roles["honest"])),
                 "af_cond_mean": float(afc.mean()), "spec_cond_mean": float(spc.mean()),
                 "strategies": {}, "per_split": []}
        for r in rows:
            block["per_split"].append({
                "af_cond": float(r["af_conditioned_auc"]),
                "spec_cond": float(r["spec_conditioned_auc"]),
                "auc": {st: {str(q): float(r["strategies"][st][str(q)]["membership_auc"])
                             for q in s61.QS} for st in ("random", "adaptive")}})
        for st in s61.STRATEGIES:
            block["strategies"][st] = {}
            for q in s61.QS:
                auc = np.array([r["strategies"][st][str(q)]["membership_auc"] for r in rows])
                t05 = np.array([r["strategies"][st][str(q)]["membership_tpr_fpr05"] for r in rows])
                g = auc - strong
                block["strategies"][st][str(q)] = {
                    "membership_auc": float(auc.mean()),
                    "membership_tpr_fpr05": float(t05.mean()),
                    "gap_vs_strongest_passive": float(g.mean()),
                    "gap_vs_strongest_passive_ci95": _boot_ci(g),
                    "pos_splits": int((g > 0).sum())}
        res["sizes"][tag] = block
        print("%s " % tag + " | ".join(
            "%s_q%d gap%+.4f[%+.4f,%+.4f]" % (
                st, q, block["strategies"][st][str(q)]["gap_vs_strongest_passive"],
                block["strategies"][st][str(q)]["gap_vs_strongest_passive_ci95"][0],
                block["strategies"][st][str(q)]["gap_vs_strongest_passive_ci95"][1])
            for st in ("random", "adaptive") for q in s61.QS), flush=True)
    # record the actual role sizes for the audit
    r0, r1 = nested_roles(groups, 0)
    res["role_sizes"] = {"per120": {k: int(v.size) for k, v in r0.items()},
                         "per180": {k: int(v.size) for k, v in r1.items()},
                         "shared_partner": int(np.intersect1d(r0["partner"], r1["partner"]).size),
                         "shared_nonmember": int(np.intersect1d(r0["nonmember"], r1["nonmember"]).size),
                         "shared_reference": int(np.intersect1d(r0["reference"], r1["reference"]).size),
                         "shared_shadow": int(np.intersect1d(r0["shadow"], r1["shadow"]).size),
                         "honest120_in_honest180": bool(
                             np.isin(r0["honest"], r1["honest"]).all())}
    with open(os.path.join(OUT, f"s62_sample_size_chr{chrom}.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s62_sample_size_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
