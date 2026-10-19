"""S64 - q-saturation cell for the selection-input spectral tomography (deep round).

The deep review (W2/R4) noted the capability question was answered only at one design
point, q <= 8; at q8 only 66-73 of the 120 honest dimensions are visible, so the
q-scaling is not shown to saturate. This script extends the SAME 5 split seeds and the
SAME machinery to q = 16, on the within-population cohorts (within:EUR and within:EAS
on chr22) and on the cross-population cohort with the headroom (super:3pop), and
reports whether the increments saturate.

The cohorts with headroom matter: at q8 the within-population attack is already at
AUC 1.000 / TPR@5% 1.000, so a flat q8->q16 gap there is partly 'cannot exceed
perfect'. super:3pop sits at AUC 0.907 / TPR@5% 0.58 at q8 and can still move.

Output: results/raw/s64_saturation_chr<C>.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s51_real_lowstruct as s51
import s61_spectral_tomography as s61
from real_multi import SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
QS_EXT = (1, 2, 4, 8, 16)


def cohorts_for(chrom, ld_window=1500, max_snps=20000):
    Ge, _, _ = s48.load_for(("EUR",), chrom, ld_window, max_snps)
    G3, _, sup3, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    G5, _, sup5, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    return {"within:EUR": (Ge, np.array(["EUR"] * len(Ge))),
            "super:3pop": (G3, sup3),
            "within:EAS": (G5[sup5 == "EAS"], sup5[sup5 == "EAS"])}


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "qs": list(QS_EXT), "cohorts": {}}
    for name, (G, groups) in cohorts_for(chrom).items():
        rows = [s61.cell(G, groups, s, qs=QS_EXT) for s in range(seeds)]

        def col(st, q, key):
            return np.array([r["strategies"][st][str(q)][key] for r in rows])

        afc = np.array([r["af_conditioned_auc"] for r in rows])
        spc = np.array([r["spec_conditioned_auc"] for r in rows])
        strong = np.maximum(afc, spc)
        blk = {"marker_count": int(G.shape[1]),
               "af_cond_mean": float(afc.mean()), "spec_cond_mean": float(spc.mean()),
               "strategies": {}, "per_split": []}
        for r in rows:
            blk["per_split"].append({
                "af_cond": float(r["af_conditioned_auc"]),
                "spec_cond": float(r["spec_conditioned_auc"]),
                "random_auc": {str(q): float(r["strategies"]["random"][str(q)]["membership_auc"])
                               for q in QS_EXT}})
        for st in s61.STRATEGIES:
            blk["strategies"][st] = {}
            for q in QS_EXT:
                auc = col(st, q, "membership_auc")
                g = auc - strong
                blk["strategies"][st][str(q)] = {
                    "membership_auc": float(auc.mean()),
                    "membership_tpr_fpr05": float(col(st, q, "membership_tpr_fpr05").mean()),
                    "gap_vs_strongest_passive": float(g.mean()),
                    "gap_vs_strongest_passive_ci95": s61._boot_ci(g),
                    "pos_splits": int((g > 0).sum()),
                    "visible_share": float(col(st, q, "visible_share").mean()),
                    "rank_visible": float(col(st, q, "rank_visible").mean()),
                    "scatter_recon_err": float(col(st, q, "scatter_recon_err").mean()),
                    "target_err_outside": float(col(st, q, "target_err_outside").mean())}
        out["cohorts"][name] = blk
        print("### %s  af_cond %.3f  spec_cond %.3f" % (name, afc.mean(), spc.mean()),
              flush=True)
        for st in s61.STRATEGIES:
            print("  %-8s " % st + " | ".join(
                "q%d auc%.3f tpr05%.3f gap%+.4f vis%.3f recon%.3f" % (
                    q, blk["strategies"][st][str(q)]["membership_auc"],
                    blk["strategies"][st][str(q)]["membership_tpr_fpr05"],
                    blk["strategies"][st][str(q)]["gap_vs_strongest_passive"],
                    blk["strategies"][st][str(q)]["visible_share"],
                    blk["strategies"][st][str(q)]["scatter_recon_err"]) for q in QS_EXT),
                flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s64_saturation_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
