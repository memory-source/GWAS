"""Synthetic algebra checks for research hypotheses, not a real genomic attack.
Run: /opt/miniforge/bin/python research-plan/20261016_deep_search/algebra_witnesses.py
Only generated arrays are used. Outputs are aggregate errors and ranks.
"""
from pathlib import Path
import json
import numpy as np

def spectral_witness():
    rng = np.random.default_rng(20261016)
    p, nh, na, k, rounds = 24, 80, 20, 3, 12
    freq = rng.uniform(.2, .5, p)
    honest = rng.binomial(2, freq, size=(nh, p))
    own_initial = rng.binomial(2, freq, size=(na, p))
    mu = np.concatenate([honest, own_initial]).mean(axis=0)
    scale = np.sqrt(mu * (1 - mu / 2))
    honest_z = (honest - mu) / scale
    honest_scatter = honest_z.T @ honest_z
    vectors, measurements, exports, checkpoints = [], [], [], []
    max_identity_error = 0.
    for t in range(rounds):
        own = own_initial.copy()
        if t:
            pair = rng.choice(na, 2, replace=False)
            block = int(rng.integers(0, 6))
            columns = np.arange(block * 4, (block + 1) * 4)
            own[np.ix_(pair, columns)] = own[np.ix_(pair[::-1], columns)]
        assert np.array_equal(np.sort(own, axis=0), np.sort(own_initial, axis=0))
        own_z = (own - mu) / scale
        known_scatter = own_z.T @ own_z
        values, basis = np.linalg.eigh(honest_scatter + known_scatter)
        v, lam = basis[:, -k:], values[-k:]
        b = v * lam - known_scatter @ v
        max_identity_error = max(max_identity_error, float(np.max(np.abs(b - honest_scatter @ v))))
        vectors.append(v)
        measurements.append(b)
        exports.append((v, lam, known_scatter))
        q = np.concatenate(vectors, axis=1)
        y = np.concatenate(measurements, axis=1)
        if t + 1 in (1, 2, 4, 8, 12):
            singular = np.linalg.svd(q, compute_uv=False)
            estimate = y @ np.linalg.pinv(q)
            checkpoints.append({
                "releases": t + 1,
                "measurement_span_rank": int(np.linalg.matrix_rank(q)),
                "smallest_reported_singular_value": float(singular[-1]),
                "relative_scatter_error": float(np.linalg.norm(estimate - honest_scatter) / np.linalg.norm(honest_scatter))
            })
    quantization = []
    for digits in (4, 6, 8):
        qs, ys = [], []
        for v, lam, known_scatter in exports:
            vq, lq = np.round(v, digits), np.round(lam, digits)
            qs.append(vq)
            ys.append(vq * lq - known_scatter @ vq)
        q, y = np.concatenate(qs, axis=1), np.concatenate(ys, axis=1)
        quantization.append({
            "decimal_places": digits,
            "releases": rounds,
            "relative_scatter_error": float(np.linalg.norm(y @ np.linalg.pinv(q) - honest_scatter) / np.linalg.norm(honest_scatter))
        })
    assert checkpoints[-1]["relative_scatter_error"] < 1e-9
    return {
        "scope": "Artificial independent-locus dosage arrays, no real people, no LD model, no membership classifier.",
        "dimensions": {"snps": p, "honest_rows": nh, "own_rows": na, "top_k": k},
        "input_change": "Each release starts from own_initial; after the baseline, exactly one 4-SNP block is exchanged between two own rows.",
        "marginal_histograms_preserved": True,
        "changed_rows_per_nonbaseline_release": 2,
        "changed_fraction_of_own_rows": 0.1,
        "max_eigen_equation_residual": max_identity_error,
        "checkpoints": checkpoints,
        "rounded_export_checks": quantization,
        "not_established": ["real genomic feasibility", "an adaptive selection algorithm", "small-query membership improvement", "real export normalization or missingness", "population or individual reconstruction"]
    }

def pair_sum_witness():
    rng = np.random.default_rng(20261017)
    n, traits = 6, 4
    residuals = rng.normal(size=(n, traits))
    residuals -= residuals.mean(axis=0)
    cycle = [(i, (i + 1) % n) for i in range(n)]
    outputs = []
    for name, edges in (("even_cycle", cycle), ("even_cycle_with_odd_chord", cycle + [(0, 2)])):
        incidence = np.zeros((len(edges), n))
        for j, (a, b) in enumerate(edges):
            incidence[j, a] = incidence[j, b] = 1.
        sums = incidence @ residuals
        estimate = np.linalg.pinv(incidence) @ sums
        outputs.append({
            "graph": name, "rows_per_measurement": 2,
            "incidence_rank": int(np.linalg.matrix_rank(incidence)),
            "individuals": n, "traits": traits,
            "max_recovery_error": float(np.max(np.abs(estimate - residuals)))
        })
    assert outputs[0]["incidence_rank"] == n - 1
    assert outputs[1]["incidence_rank"] == n
    assert outputs[1]["max_recovery_error"] < 1e-12
    return {
        "scope": "A known-incidence linear inverse subproblem, not blind carrier discovery or a real GWAS release.",
        "results": outputs,
        "not_established": ["recovering an unknown carrier graph", "recovering raw clinical phenotypes", "sufficient overlap in real local genealogies", "handling GLMM weighting, missingness, or MAC filters"]
    }

def main():
    result = {
        "date": "2026-10-16", "seeded_and_synthetic": True,
        "spectral": spectral_witness(), "pair_sums": pair_sum_witness()
    }
    destination = Path(__file__).with_name("algebra_witnesses.json")
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
