"""1000G Phase-3 IMPUTE2 genetic map (bp -> centiMorgan).

The map files are whitespace-delimited numeric tables (three columns:
position, COMBINED_rate(cM/Mb), Genetic_Map(cM)) downloaded to
data/external/genetic_map/genetic_map_chr<C>.txt. `numpy.loadtxt` parses them;
`bp_to_cm` linearly interpolates the cumulative cM at an arbitrary bp, so a
marker's physical position on chr1/chr2/chr22 can be expressed in genetic
distance and an exact-match run can be measured in cM.
"""
import os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAPDIR = os.path.join(ROOT, "data", "external", "genetic_map")


def map_path(chrom):
    return os.path.join(MAPDIR, f"genetic_map_chr{chrom}.txt")


def load_map(chrom):
    """Return (pos_bp (M,), cM (M,)) sorted by position."""
    p = map_path(chrom)
    if not os.path.exists(p):
        raise FileNotFoundError(f"genetic map missing: {p}")
    tab = np.loadtxt(p, skiprows=1)
    order = np.argsort(tab[:, 0])
    return tab[order, 0], tab[order, 2]


def bp_to_cm(pos_query, chrom):
    """Interpolate cumulative cM at each queried bp position."""
    pos, cm = load_map(chrom)
    return np.interp(np.asarray(pos_query, dtype=np.float64), pos, cm)


if __name__ == "__main__":
    for c in (1, 2, 22):
        pos, cm = load_map(c)
        print(f"chr{c}: {pos.size} rows, {pos[0]/1e6:.2f}-{pos[-1]/1e6:.2f} Mb, "
              f"{cm[-1]:.1f} cM")
