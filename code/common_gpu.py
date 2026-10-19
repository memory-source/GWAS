"""common_gpu.py - torch port of the leakage-measurement kernels (GPU).

Mirrors code/common.py so results are comparable. Runs under the SYSTEM python3
(torch 2.2.0+cu121), not the project conda env.

Memory model (why panels are HOST-resident):
  A biobank-scale panel is n x M int8, e.g. 9e4 x 1e6 = 90 GB: it cannot live on
  a 24 GB GPU. So the simulator produces a HOST numpy int8 array; the PCA (Gram)
  and the attack scores STREAM marker blocks from host to GPU. The peak GPU
  footprint is then O(n^2) (the Gram) or O(n*blk), independent of M.

Precision: TF32 is left OFF -- the PC channel's magnitude is defined by the
*fitted* loadings, so TF32 would corrupt them.
"""
import numpy as np
import torch

DEV = torch.device("cuda" if torch.cuda.is_available() else "cpu")
torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False


# ---------------------------------------------------------------- simulator ---
def simulate_population_host(n_groups=3, n_per_group=1000, n_snps=4000, fst=0.08,
                             ld_blocks=100, ld_size=6, seed=0, blk=None, budget=2e9,
                             pct="beta"):
    """Balding-Nichols K-population panel + adjacent LD injection, on HOST.

    Returns dict(G=np.int8 (n, M), ancestry=np.int64 (n,), n_groups). Memory
    O(n*M) bytes for G, and O(n*blk) for the transient float32 success probs --
    marker columns are drawn in blocks so a biobank panel (n*M up to ~1e11) never
    materialises a full float32 intermediate (which would be 4x the int8 panel).
    ``blk`` defaults to a block size that keeps that transient under ``budget``
    bytes regardless of n.

    ``pct`` selects ONLY the per-population allele-frequency draw:
      'beta'  -- Balding-Nichols gf ~ Beta(a, b). Exact, but the Beta draw is the
                 dominant cost at large cells (~5-6x two random() draws).
      'gauss' -- diffusion-limit approximation gf = p + N(0, sqrt(p(1-p) Fst)),
                 clipped to [0,1]; matches BN to O(Fst^2). Used for the
                 biobank-scale grid where the Beta draw, not the genotype
                 sampling, is the bottleneck."""
    rng = np.random.default_rng(seed)
    n = n_groups * n_per_group
    if blk is None:
        blk = int(max(2048, min(200000, budget / (4.0 * n))))
    anc = np.repeat(np.arange(n_groups), n_per_group)
    base = np.clip(rng.uniform(0.05, 0.95, size=n_snps), 1e-4, 1 - 1e-4)
    fst = max(fst, 1e-6)
    a = base * (1 - fst) / fst
    b = (1 - base) * (1 - fst) / fst
    sd = np.sqrt(base * (1 - base) * fst).astype(np.float32)
    G = np.empty((n, n_snps), dtype=np.int8)
    for s in range(0, n_snps, blk):
        e = min(s + blk, n_snps)
        if pct == "gauss":
            gf = base[None, s:e] + rng.normal(size=(n_groups, e - s)) * sd[None, s:e]
            gf = np.clip(gf, 1e-6, 1 - 1e-6).astype(np.float32)           # (K, blk)
        else:
            gf = rng.beta(a[None, s:e], b[None, s:e], size=(n_groups, e - s))
        F = gf[anc]                                                       # (n, blk)
        # Binomial(2, p) via two Bernoulli draws (much faster than the
        # Generator.binomial slow path). NB: cast each draw to int8 before adding
        # -- numpy `+` on bool arrays is logical-OR.
        g1 = (rng.random(F.shape, dtype=np.float32) < F).astype(np.int8)
        G[:, s:e] = g1 + (rng.random(F.shape, dtype=np.float32) < F).astype(np.int8)
        del F, g1, gf
    if ld_blocks and ld_size:
        starts = rng.choice(max(1, n_snps - ld_size), size=min(ld_blocks, n_snps),
                            replace=False)
        latent_hi = rng.random(n) > 0.5
        for s in starts:
            cols = slice(int(s), min(int(s) + ld_size, n_snps))
            flip = (rng.random(n) < 0.15) & latent_hi
            G[flip, cols] = np.clip(G[flip, cols].astype(np.int16) + 1, None, 2).astype(np.int8)
    return dict(G=G, ancestry=anc, n_groups=n_groups)


def split_population3_host(pop, n_site=300, n_ref=300, n_non=300, seed=0):
    """Per-group 3-way disjoint split (matches CPU split_population3)."""
    G, anc = pop["G"], pop["ancestry"]
    K = pop["n_groups"]; per = G.shape[0] // K
    rng = np.random.default_rng(seed)
    s_i, r_i, n_i = [], [], []
    for k in range(K):
        idx = rng.permutation(np.arange(k * per, (k + 1) * per))
        a = min(n_site, per // 3); b = min(n_ref, per // 3); c = min(n_non, per // 3)
        s_i.append(idx[:a]); r_i.append(idx[a:a + b]); n_i.append(idx[a + b:a + b + c])
    take = lambda ii: (np.concatenate(ii))
    si, ri, ni = take(s_i), take(r_i), take(n_i)
    mk = lambda ii: dict(G=G[ii], ancestry=anc[ii], n_groups=K)
    return mk(si), mk(ri), mk(ni)


# -------------------------------------------------------------------- PCA -----
def gram_svd_host(G_host, indices=None, k=10, block=20000, device=DEV):
    """Exact top-k PCA of a centered host int8 panel via the n x n Gram matrix.

    Streams marker blocks host->GPU, so the GPU holds only the (n x n) Gram.
    Returns (eig (k,) desc on host, loadings (k, M) on GPU)."""
    G = G_host if indices is None else G_host[indices]
    n, M = G.shape
    mu = G.mean(0, dtype=np.float64)                     # (M,) host
    Gm = torch.zeros(n, n, device=device)
    for s in range(0, M, block):
        e = min(s + block, M)
        B = torch.as_tensor(G[:, s:e], device=device, dtype=torch.float32)
        B = B - torch.as_tensor(mu[s:e], device=device, dtype=torch.float32)[None]
        Gm.add_(B @ B.T)
        del B
    ev, U = torch.linalg.eigh(Gm)                        # ascending
    idx = torch.argsort(ev, descending=True)[:k]
    ev_k = ev[idx].clamp(min=1e-9)
    U_k = U[:, idx]
    S_k = torch.sqrt(ev_k)
    L = torch.empty(k, M, device=device)
    for s in range(0, M, block):
        e = min(s + block, M)
        B = torch.as_tensor(G[:, s:e], device=device, dtype=torch.float32)
        B = B - torch.as_tensor(mu[s:e], device=device, dtype=torch.float32)[None]
        L[:, s:e] = (B.T @ (U_k / S_k[None, :])).T
        del B
    del Gm, U
    torch.cuda.empty_cache()
    return ev_k.detach().cpu().numpy(), L.detach()


# ----------------------------------------------------------------- attacks ----
def _log_tables(f):
    f = f.clamp(1e-6, 1 - 1e-6)
    return torch.stack([torch.log((1 - f) ** 2),
                        torch.log(2 * f * (1 - f)),
                        torch.log(f * f)], 0)


def af_lrt_host(G_rows_host, af_s_host, af_r_host, block=50000, device=DEV):
    """AF-LRT scores streamed per marker block. af_s/af_r (M,) host. -> (n,) GPU."""
    n = G_rows_host.shape[0]; M = G_rows_host.shape[1]
    out = torch.zeros(n, device=device)
    AFS = torch.as_tensor(af_s_host, device=device).float()
    AFR = torch.as_tensor(af_r_host, device=device).float()
    for s in range(0, M, block):
        e = min(s + block, M)
        g = torch.as_tensor(G_rows_host[:, s:e], device=device).long()
        Ls = _log_tables(AFS[s:e][None])
        Lr = _log_tables(AFR[s:e][None])
        d = (Ls - Lr)[:, 0]                             # (3, blk)
        for c in range(3):
            out += d[c] @ (g == c).float().t()
        del g, d, Ls, Lr
    return out


def af_lrt_gpu(G_gpu, af_s, af_r, mask=None):
    """In-GPU AF-LRT for small panels (V-many af tables supported)."""
    g = G_gpu.long()
    if mask is not None:
        g = g[:, mask]; af_s = af_s[..., mask]; af_r = af_r[..., mask]
    single = af_s.dim() == 1
    if single:
        af_s = af_s[None, :]; af_r = af_r[None, :]
    Ls, Lr = _log_tables(af_s.float()), _log_tables(af_r.float())
    V, n = af_s.shape[0], g.shape[0]
    out = torch.zeros(V, n, device=g.device)
    for c in range(3):
        out += (Ls[c] - Lr[c]) @ (g == c).float().t()
    return out[0] if single else out


def proj_dist_host(G_rows_host, loadings_gpu, mean_host, block=50000, device=DEV):
    """Projection distance to the study centroid, streamed. loadings (k, M) GPU,
    mean (M,) host. -> (n,) GPU."""
    n = G_rows_host.shape[0]; M = G_rows_host.shape[1]
    k = loadings_gpu.shape[0]
    proj = torch.zeros(n, k, device=device)
    for s in range(0, M, block):
        e = min(s + block, M)
        B = torch.as_tensor(G_rows_host[:, s:e], device=device, dtype=torch.float32)
        B = B - torch.as_tensor(mean_host[s:e], device=device, dtype=torch.float32)[None]
        proj += B @ loadings_gpu[:, s:e].t()
        del B
    return proj.norm(dim=1) / (M ** 0.5)


def proj_dist_gpu(G_gpu, loadings, mean, mask=None, k=None):
    if mask is not None:
        G_gpu = G_gpu[:, mask]; mean = mean[..., mask]
        loadings = loadings[..., mask] if loadings.dim() == 3 else loadings[:, mask]
    X = G_gpu.to(torch.float32) - mean.to(torch.float32)[None]
    if loadings.dim() == 3:
        L = loadings if k is None else loadings[:, :k]
        proj = torch.einsum('nm,vkm->vkn', X, L); d = proj.norm(dim=1)
    else:
        L = loadings if k is None else loadings[:k]
        d = (X @ L.t()).norm(dim=1)
    return d / (X.shape[1] ** 0.5)


# ------------------------------------------------------------------ metrics ---
def auroc_gpu(scores, labels):
    s = scores.flatten(); y = labels.flatten().to(device=s.device, dtype=torch.float32)
    order = torch.argsort(s)
    ranks = torch.empty_like(s)
    ranks[order] = torch.arange(1, s.numel() + 1, device=s.device, dtype=torch.float32)
    npos = y.sum(); nneg = s.numel() - npos
    return ((ranks * y).sum() - npos * (npos + 1) / 2) / (npos * nneg)


def batch_auroc_gpu(scores, label_batch):
    s = scores.flatten()
    order = torch.argsort(s)
    ranks = torch.empty_like(s)
    ranks[order] = torch.arange(1, s.numel() + 1, device=s.device, dtype=torch.float32)
    y = label_batch.to(device=s.device, dtype=torch.float32)
    sr = (ranks[None, :] * y).sum(1)
    npos = y.sum(1); nneg = s.numel() - npos
    return (sr - npos * (npos + 1) / 2) / (npos * nneg)


def stacked_split(site_G, non_G):
    X = torch.cat([site_G, non_G], 0)
    y = torch.cat([torch.ones(site_G.shape[0], device=X.device),
                   torch.zeros(non_G.shape[0], device=X.device)])
    return X, y


# ---- backwards-compatible in-GPU simulator (small panels only) ----
def simulate_population_gpu(n_groups=3, n_per_group=1000, n_snps=4000, fst=0.08,
                            ld_blocks=100, ld_size=6, seed=0, device=DEV):
    p = simulate_population_host(n_groups, n_per_group, n_snps, fst, ld_blocks, ld_size, seed)
    return dict(G=torch.as_tensor(p["G"], device=device),
                ancestry=torch.as_tensor(p["ancestry"], device=device),
                n_groups=n_groups)


def split_population3_gpu(pop, n_site=300, n_ref=300, n_non=300, seed=0):
    G, anc = pop["G"], pop["ancestry"]; K = pop["n_groups"]; per = G.shape[0] // K
    rng = np.random.default_rng(seed)
    s_i, r_i, nm_i = [], [], []
    for k in range(K):
        idx = rng.permutation(np.arange(k * per, (k + 1) * per))
        a = min(n_site, per // 3); b = min(n_ref, per // 3); c = min(n_non, per // 3)
        s_i.append(idx[:a]); r_i.append(idx[a:a + b]); nm_i.append(idx[a + b:a + b + c])
    to = lambda ii: G[torch.as_tensor(np.concatenate(ii), device=G.device)]
    ta = lambda ii: anc[torch.as_tensor(np.concatenate(ii), device=G.device)]
    return (dict(G=to(s_i), ancestry=ta(s_i), n_groups=K),
            dict(G=to(r_i), ancestry=ta(r_i), n_groups=K),
            dict(G=to(nm_i), ancestry=ta(nm_i), n_groups=K))
