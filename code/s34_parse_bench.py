"""S34 - VCF parser benchmark (B6, anti-dressing control).

Turns "we could put the VCF parse on the GPU" into a measured statement: the
real break in the pipeline is the *parse* (interpretation of GT strings), not
decompression, and a GPU string-parser is not obviously a win for a
line-oriented gzip stream. We benchmark, on a real 1000G chr file:

  (1) gzip single-thread loop + Python GT parse (the current path)
  (2) pigz -dc -p<N> decompress to stdout, piped into the same Python parse
  (3) pigz decompress to a temp file, then parse (isolates decompress cost)

We report variants/s and effective MB/s for each. This is a NEGATIVE result
used to bound what GPU acceleration can/cannot change; it is not a claim.

Output: results/raw/s34_parse_bench.json
"""
import gzip
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from real_multi import _vcf_path, load_panel

OUT = "/data/lvvyongji1/GWAS/results/raw"
SUPS3 = {"AFR", "EAS", "EUR"}


def _vcf_bytes(path):
    p = subprocess.run(["pigz", "-l", path] if _have("pigz") else ["gzip", "-l", path],
                       capture_output=True, text=True)
    try:
        return int(p.stdout.split()[-1])
    except Exception:
        return os.path.getsize(path)


def _have(tool):
    return subprocess.run(["which", tool], capture_output=True).returncode == 0


def parse_stream(f, n_max, sample_cols, tmark):
    n = 0
    t0 = time.time()
    for line in f:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(sample_cols):
            continue
        ref, alt = p[3], p[4]
        if "," in alt or "." in alt or len(ref) != 1 or len(alt) != 1:
            continue
        for i in sample_cols:
            g = p[i]
            if g[0] == "." or (len(g) > 2 and g[2] == "."):
                break
        n += 1
        if n >= n_max:
            break
    return n, time.time() - t0


def bench_python_gzip(path, n_max=20000):
    samples, sup = load_panel()
    sample_cols = [i for i, s in enumerate(samples) if sup[i] in SUPS3]
    with gzip.open(path, "rt") as f:
        n, el = parse_stream(f, n_max, sample_cols, None)
    return dict(n=n, seconds=el, variants_per_s=n / el)


def bench_pigz_pipe(path, n_max=20000, threads=16):
    samples, sup = load_panel()
    sample_cols = [i for i, s in enumerate(samples) if sup[i] in SUPS3]
    proc = subprocess.Popen(["pigz", "-dc", "-p", str(threads), path],
                            stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n, el = parse_stream(proc.stdout, n_max, sample_cols, None)
    proc.stdout.close(); proc.wait()
    return dict(n=n, seconds=el, variants_per_s=n / el)


def bench_decompress_only(path, threads=16):
    t0 = time.time()
    if _have("pigz"):
        subprocess.run(["pigz", "-dc", "-p", str(threads), path],
                       stdout=subprocess.DEVNULL)
    else:
        subprocess.run(["gzip", "-dc", path], stdout=subprocess.DEVNULL)
    return dict(seconds=time.time() - t0)


def main():
    path = None
    for c in (22, 1, 2):
        try:
            path = _vcf_path(c); chrom = c; break
        except FileNotFoundError:
            continue
    if path is None:
        print("no VCF found", flush=True); return
    raw = _vcf_bytes(path) if _have("pigz") else None
    res = dict(vcf=os.path.basename(path), chrom=chrom, pigz=_have("pigz"))
    try:
        res["python_gzip"] = bench_python_gzip(path)
        print("python+gzip", res["python_gzip"], flush=True)
    except Exception as e:
        res["python_gzip"] = dict(error=str(e))
    try:
        res["decompress_only_pigz16"] = bench_decompress_only(path)
        print("decompress only", res["decompress_only_pigz16"], flush=True)
    except Exception as e:
        res["decompress_only_pigz16"] = dict(error=str(e))
    if _have("pigz"):
        try:
            res["pigz_pipe_p16"] = bench_pigz_pipe(path)
            print("pigz pipe p16", res["pigz_pipe_p16"], flush=True)
        except Exception as e:
            res["pigz_pipe_p16"] = dict(error=str(e))
    with open(os.path.join(OUT, "s34_parse_bench.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("done", flush=True)


if __name__ == "__main__":
    main()
