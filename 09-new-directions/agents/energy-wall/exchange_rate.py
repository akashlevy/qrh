"""
Exponent bookkeeping for the reflected marked energy beyond the wall (P78 Lemma hybrid-energy,
binding dyad O = 0, l_b = 0, v = T_d).  All quantities are log_Z.

  rows    K = Z^M          (residual squarefree rows k, quadratic characters psi_k = chi_k^3)
  columns N = Z^(2M+l-1)   (reflected squarefree dual index n, coefficient beta(n) X(n))
  marks   L = Z^l          (product of active slot primes; dual twist chi_P(n)^{-2})
  wall excess  w = M + l - 1 = log_Z(N/K)  (energy exponent minus diagonal M)

"Transfer" = Cauchy-Schwarz a coherent variable V (log-length v) out of the column sum into the
row family, then assume the OPTIMAL (diagonal, no (CN)^{2/3} term) large sieve for the enlarged
family.  The resulting excess over the diagonal is computed for each kind of variable.
Conclusion (for every v): no transfer beats w whenever l >= w (equivalently M <= 1).
"""
import numpy as np

def excesses(M, l):
    w = M + l - 1
    v = np.linspace(0, 1.0, 2001)
    out = {}
    out["paper (quadratic LS column term)"] = w
    # slot of the marks moved to rows (needs KP_1 >= N for diagonality): cost v, columns unchanged
    vs = v[v <= l]
    out["slot transfer, min over v<=l"] = float(np.min(np.maximum(0, w - vs) + vs))
    # theta factor n2 of the dual index (twisted multiplicativity), coupled to the marks:
    # columns become m = n1 P^2 (dense count N L^2 / N2) or sparse count N L / N2 (optimistic)
    out["TM factor (mark-coupled, dense sieve)"] = float(np.min(np.maximum(v, w + 2 * l - v)))
    out["TM factor (mark-coupled, optimistic sparse sieve)"] = float(np.min(np.maximum(v, w + l - v)))
    # theta factor AND marks both moved to rows
    out["TM factor + marks to rows"] = float(np.min(np.maximum(0, w - 2 * v - l) + v + l))
    # hypothetical decoupled factor of the dual index (does not exist at the binding dyad)
    out["decoupled factor (hypothetical)"] = float(np.min(np.maximum(0, w - 2 * v) + v))
    # unreflected (primal) side: sextic LS, conductor range Z^{1+l}, BGL term (K N)^{2/3}
    out["primal, BGL sextic sieve"] = max(1 + l - M, (2 + 2 * l - M) / 3 - 0)
    out["primal, if sparsity (Z effective columns) were exploitable"] = 1 - M
    return w, out

if __name__ == "__main__":
    for (M, l, tag) in [(0.36728992632054164 + 0.49503135088645195, 0.20767872279300653, "model optimum at wE=0.07"),
                        (5 / 6, 1 / 6, "P78 geometry (on the wall)")]:
        w, out = excesses(M, l)
        print(f"{tag}: M={M:.4f} l={l:.4f} w=M+l-1={w:.4f}")
        for k, x in out.items():
            print(f"   {k:58s} excess {x:+.4f}   {'(<= w: gain)' if x < w - 1e-9 else ''}")
