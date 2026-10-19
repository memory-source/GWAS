"""Numerically verify fixed-subspace moment identities on synthetic genotypes.
This is an algebra check, not a membership-attack experiment or novelty proof.
"""
from pathlib import Path
import hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np

def residual_moments(n, V, eigenvalues, known_raw_G, released_mean, released_scale):
    """Attacker view: receives only release parameters and known member rows."""
    u=n-len(known_raw_G)
    z_known=((known_raw_G-released_mean)/released_scale)@V
    remaining_mean=-z_known.sum(0)/u
    remaining_scatter=n*np.diag(eigenvalues)-z_known.T@z_known-u*np.outer(remaining_mean,remaining_mean)
    return remaining_mean,remaining_scatter

def main():
    rng=np.random.default_rng(20261015)
    n,p,k=48,24,5
    f=rng.uniform(0.1,0.9,p)
    G=rng.binomial(2,f,size=(n,p)).astype(float)
    mu=G.mean(0)
    af=G.mean(0)/2
    scale=np.sqrt(2*af*(1-af))
    X=(G-mu)/scale
    ev,V=np.linalg.eigh(X.T@X/n)
    order=np.argsort(ev)[::-1][:k]
    ev,V=ev[order],V[:,order]
    order_rows=rng.permutation(n)
    cells=[]
    for m in (0,12,36,44):
        known=order_rows[:m]
        unknown=order_rows[m:]
        inferred_mean,inferred_scatter=residual_moments(n,V,ev,G[known],mu,scale)
        # The evaluator alone reads unknown member rows to check the identity.
        z_unknown=X[unknown]@V
        direct_mean=z_unknown.mean(0)
        centered=z_unknown-direct_mean
        direct_scatter=centered.T@centered
        cells.append(dict(n=n,p=p,k=k,known_members=m,unknown_members=n-m,
                          max_abs_mean_error=float(np.max(np.abs(inferred_mean-direct_mean))),
                          max_abs_scatter_error=float(np.max(np.abs(inferred_scatter-direct_scatter))),
                          relative_scatter_frobenius_error=float(np.linalg.norm(inferred_scatter-direct_scatter)/max(np.linalg.norm(direct_scatter),1e-15)),
                          inferred_scatter_smallest_eigenvalue=float(np.linalg.eigvalsh(inferred_scatter).min()),
                          direct_scatter_rank=int(np.linalg.matrix_rank(direct_scatter,tol=1e-9)),
                          rank_upper_bound=min(k,n-m-1)))
    Q,_=np.linalg.qr(rng.normal(size=(k,k)))
    norms=np.linalg.norm(X@V,axis=1)
    rotated=np.linalg.norm(X@V@Q,axis=1)
    output=dict(scope="Synthetic exact-arithmetic identity up to floating-point error; no real payload, attack AUC, privacy guarantee, or new-theorem claim.",
                seed=20261015,normalization="S=X.T@X/n; X=(G-released_mean)/released_scale",
                input_genotype_sha256=hashlib.sha256(G.tobytes()).hexdigest(),
                mean_equals_twice_AF_max_error=float(np.max(np.abs(mu-2*af))),
                rotation_projection_norm_max_error=float(np.max(np.abs(norms-rotated))),cells=cells)
    assert max(c["max_abs_scatter_error"] for c in cells)<1e-9
    assert max(c["max_abs_mean_error"] for c in cells)<1e-10
    assert output["rotation_projection_norm_max_error"]<1e-10
    path=Path(__file__).with_name("algebra_check.json")
    path.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps(output,indent=2))
if __name__=="__main__":
    main()
