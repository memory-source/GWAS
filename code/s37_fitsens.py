"""S37 wrap-up - PCA fit-size sensitivity sweep at the largest cells.

The s37 cell fits the study-fitted PC basis on at most `n_fit` study rows (the
GPU Gram is n_fit x n_fit). At nsp=8000/16000 the study has 24000/48000 rows, so
a cell reported with n_fit=20000 uses a SUBSAMPLED fit. This probe sweeps the fit
size at the two largest cells to check whether the AF-first conclusion depends on
it:

  nsp=8000  (n_site=24000): n_fit in {5000, 10000, 20000, full}
  nsp=16000 (n_site=48000): n_fit in {5000, 10000, 20000}

`full` (n_fit=n_site) is run for the nsp=8000 cell (n_site=24000); for the
nsp=16000 cell (n_site=48000) the full fit is not run. The gap is compared to the
per-cell std.

Output: results/raw/s37_fitsens.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
import s37_realsynth_scale as s37

OUT = "/data/lvvyongji1/GWAS/results/raw"
SEEDS = (0, 1, 2, 3, 4)
SWEEP = {
    (8000, 8000): (5000, 10000, 20000, 24000),
    (16000, 3000): (5000, 10000, 20000),
}


def run():
    out = {}
    for (nsp, M), fits in SWEEP.items():
        key = f"nsp{nsp}_m{M}"
        cell = {}
        for nf in fits:
            rows = [s37.cell(nsp, M, s, n_fit=nf) for s in SEEDS]
            cell[f"nfit{nf}"] = dict(
                n_fit=int(min(nf, rows[0]["n_site"])),
                fit_frac=float(rows[0]["fit_frac"]),
                af_only=float(np.mean([r["af_only"] for r in rows])),
                pc_calibrated=float(np.mean([r["pc_calibrated"] for r in rows])),
                delta_pc_minus_af=float(np.mean([r["delta_pc_minus_af"] for r in rows])),
                delta_pc_minus_af_std=float(np.std([r["delta_pc_minus_af"] for r in rows])),
                n_site=int(rows[0]["n_site"]), n_seeds=len(rows))
            c = cell[f"nfit{nf}"]
            print("%s nfit %5d  fit_frac %.3f  AF %.3f  PC %.3f  gap %+.3f (sd %.3f)" % (
                key, c["n_fit"], c["fit_frac"], c["af_only"], c["pc_calibrated"],
                c["delta_pc_minus_af"], c["delta_pc_minus_af_std"]), flush=True)
            out[key] = cell
            with open(os.path.join(OUT, "s37_fitsens.json"), "w") as f:
                json.dump(out, f, indent=2)
    print("fitsens sweep done", flush=True)
    return out


if __name__ == "__main__":
    run()
