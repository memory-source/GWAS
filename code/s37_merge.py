"""Merge the round-11 sharded results into single files + a consolidated table.

s37 (biobank semi-synthetic) wrote one JSON per shard; this pools the shards into
results/raw/s37_realsynth_scale.json and prints the boundary table alongside the
s27 parametric grid so the paper can cite one place.

Output: results/raw/s37_realsynth_scale.json (merged),
        refine-logs/ROUND11_RESULTS.md (table).
"""
import glob
import json
import os

RAW = "/data/lvvyongji1/GWAS/results/raw"
LOGS = "/data/lvvyongji1/GWAS/refine-logs"


def merge_s37():
    out = dict(grid={}, controls={})
    for f in sorted(glob.glob(os.path.join(RAW, "s37_realsynth_scale_shard*.json"))):
        d = json.load(open(f))
        for k, v in d.get("grid", {}).items():
            out["grid"][k] = v
        for k, v in d.get("controls", {}).items():
            out["controls"][k] = v
    out["n_shards"] = len(glob.glob(os.path.join(RAW, "s37_realsynth_scale_shard*.json")))
    with open(os.path.join(RAW, "s37_realsynth_scale.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


def boundary_table(merged):
    keys = ["af_only", "pc_calibrated", "combined", "delta_pc_minus_af",
            "pc_refbasis", "pc_eigonly", "pc_in_sample_gap", "fst_like",
            "ev_top10_frac"]
    rows = []
    order = sorted(merged["grid"], key=lambda k: int(k.split("_")[0][3:]))
    for k in order:
        g = merged["grid"][k]
        rows.append((k, {kk: (g.get(kk) if kk in g else float("nan")) for kk in keys},
                     g.get("n_seeds_used")))
    return rows


def main():
    m = merge_s37()
    rows = boundary_table(m)
    lines = ["# Round 11 — semi-synthetic biobank cohort (s37)",
             "",
             "Cohorts mosaic-copied (Li–Stephens) from the real 1000G phased library; "
             "same attack battery as the parametric s27 grid.",
             "",
             "| cell | af_only | pc_calibrated | delta(pc-af) | pc_refbasis | pc_eigonly | "
             "in_sample_gap | fst_like | n_seeds |",
             "|------|---------|---------------|--------------|-------------|------------|"
             "----------------|----------|---------|"]
    for k, r, ns in rows:
        lines.append("| %s | %.3f | %.3f | %+.3f | %.3f | %.3f | %.3f | %.3f | %s |" % (
            k, r["af_only"], r["pc_calibrated"], r["delta_pc_minus_af"],
            r["pc_refbasis"], r["pc_eigonly"], r["pc_in_sample_gap"], r["fst_like"], ns))
    lines.append("")
    lines.append("## Controls (LD decay bp-bins | relatedness | eig10 frac)")
    lines.append("")
    lines.append("| cell | LD cohort | LD donor | relatedness cohort | donor | eig10 c | d |")
    lines.append("|------|-----------|----------|--------------------|-------|---------|---|")
    for k in sorted(m["controls"], key=lambda k: int(k.split("_")[0][3:])):
        c = m["controls"][k]
        lines.append("| %s | %s | %s | %.1f | %.1f | %.3f | %.3f |" % (
            k,
            "[" + ", ".join("%.3f" % x for x in c.get("ld_cohort", [])) + "]",
            "[" + ", ".join("%.3f" % x for x in c.get("ld_donor", [])) + "]",
            c.get("relatedness_cohort", [float("nan")])[0],
            c.get("relatedness_donor", [float("nan")])[0],
            c.get("eig10_frac_cohort", float("nan")),
            c.get("eig10_frac_donor", float("nan"))))
    with open(os.path.join(LOGS, "ROUND11_RESULTS.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return m


if __name__ == "__main__":
    main()
