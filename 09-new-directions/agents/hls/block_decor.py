"""
hls: block structure of the hybrid Gram matrix (why every block-wise argument loses).

G (KL x KL) has K x K blocks G^{PP'} = Psi diag(eta_P conj eta_P') Psi^T.
  diagonal blocks  : quadratic large sieve (Heath-Brown), ||G^{PP}|| <= K + N (PROVED);
  off-diagonal     : twisted quadratic Gram, ||G^{PP'}|| <= min(K + N, KL + N/L) (PROVED: Cauchy-Schwarz / Poisson);
  block-Schur      : ||G|| <= max_P sum_P' ||G^{PP'}||  (what any block-wise proof gives);
  truth            : ||G|| ~ KL + N needs coherent cancellation ACROSS blocks.
Prints: block norms vs 2 sqrt(K' N') (random model), block-Schur sum vs the true ||G||, and the same for the
'prod' random control.
Usage: nice python3 block_decor.py N0 K L0 [nP] out.json
"""
import sys, json, math, time
import numpy as np
from hls_lib import Family, rand_control


def analyse(F, tag):
    t0 = time.time()
    Pm = F.Psi.astype(np.float64)
    nk, nP, nc = F.nk, F.nP, F.nc
    B = np.zeros((nP, nP))
    blocks = {}
    for i in range(nP):
        for j in range(i, nP):
            w = F.H[i] * np.conj(F.H[j])
            M = (Pm * w.real) @ Pm.T + 1j * ((Pm * w.imag) @ Pm.T)
            s = np.linalg.norm(M, 2)
            B[i, j] = B[j, i] = s
            if i == 0 and j <= 1:
                blocks[(i, j)] = M
    Gfull_top = F.top(k=1)[0][0]
    diag = np.diag(B); off = B[~np.eye(nP, dtype=bool)]
    schur = np.max(B.sum(axis=1))
    schur_off = np.max(B.sum(axis=1) - diag)
    # norm of the off-diagonal-block part G_off = G - blockdiag: Lanczos on the row side
    from scipy.sparse.linalg import LinearOperator, eigsh
    def mv(y):
        Y = y.reshape(nk, nP)
        Z = F.S(F.Sh(Y))                         # G Y
        # subtract block-diagonal part: for each P, G^{PP} Y[:,P] = Psi diag(|eta_P|^2) Psi^T Y[:,P]
        D = np.empty_like(Z)
        for p in range(nP):
            w = np.abs(F.H[p]) ** 2
            D[:, p] = Pm @ (w * (Pm.T @ Y[:, p]))
        return (Z - D).ravel()
    op = LinearOperator((nk * nP, nk * nP), matvec=mv, dtype=complex)
    ev = eigsh(op, k=2, which="LM", tol=1e-5)[0]
    off_norm = float(np.max(np.abs(ev)))
    res = dict(tag=tag, nk=nk, nP=nP, nc=nc, lam_G=float(Gfull_top),
               diag_block_mean=float(diag.mean()), diag_block_max=float(diag.max()),
               off_block_mean=float(off.mean()), off_block_max=float(off.max()),
               rand_model_offblock=float(2 * math.sqrt(nk * nc)),
               trivial_offblock=float(nk + nc), poisson_offblock=float(nk * nP + nc / nP),
               block_schur=float(schur), block_schur_off=float(schur_off), Goff_norm=off_norm,
               Goff_eigs=[float(x) for x in ev], time=time.time() - t0)
    print(f"[{tag}] k'={nk} P'={nP} cols={nc}: ||G||={Gfull_top:.0f} (rows+cols={nk*nP+nc})  "
          f"diag blocks {diag.mean():.0f} (max {diag.max():.0f}; K'+N'={nk+nc})  "
          f"off blocks mean {off.mean():.0f} max {off.max():.0f} (2sqrt(K'N')={2*math.sqrt(nk*nc):.0f})  "
          f"block-Schur sum {schur:.0f}  ||G_off||={off_norm:.0f}  [{time.time()-t0:.0f}s]", flush=True)
    return res


if __name__ == "__main__":
    N0, K, L0 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    nP = int(sys.argv[4]) if len(sys.argv) > 5 else None
    out = sys.argv[-1]
    F = Family(N0, K, L0, nP=nP, gauss=False)
    R = [analyse(F, "arith"), analyse(rand_control(F, "prod", seed=5), "prod-random")]
    json.dump(R, open(out, "w"), indent=1)
