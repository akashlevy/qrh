"""LP check for Attack B (B_bilinear.md).  Same model as wf/lp_const_saving.py (P^2/Y Gram piece dropped by LoA6,
slot-free low side of PARTII), theta piece g = b/6 - 2 s(M') with:
  'const' s0           : uniform saving s0 (Z-exponent);  s0 = 0 is what Attack B delivers.
  'zexp'               : the (invalid, optimistic) single-factor theta-expansion of zeta_m (Sec. 4 of B_bilinear.md):
                         ratio^2 = R/Q + R^2/Q^2 + R^2/Q^{5/3},  R = L Q^{1/2}, i.e.
                         -2s = max(3M'/2 - 1, 3M' - 2, 10M'/3 - 2)   (Z-exponents, Q = Z^{M'}, L = Z^{2M'-1}).
theta_true = min(theta_CS, theta_route) by distributivity (each route is an alternative bound)."""
import sys, inspect
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import numpy as np
import theta_split_model as T
import model_r4 as M
from scipy.optimize import minimize_scalar
src = inspect.getsource(T.lp_for_S)
src = src.replace('def lp_for_S(S, cfg, mode, return_geom=False, Rcache=None):',
                  'def lp2(S, cfg, mode, return_geom=False, Rcache=None, s0=0.0, dropq=True):')
src = src.replace('        gpieces.append(add(B, LY, d, coef=(2, -1, 1)))',
                  '        if not dropq: gpieces.append(add(B, LY, d, coef=(2, -1, 1)))')
src = src.replace('        if mode == "none":',
                  '        if mode == "const":\n'
                  '            gpieces.append(add(B, ONE, coef=(1 / 6, -2 * s0)))\n'
                  '        elif mode == "zexp":\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 1.5, -1)))\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 3, -2)))\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 10 / 3, -2)))\n'
                  '        elif mode == "none":')
ns = dict(T.__dict__); exec(src, ns); lp2 = ns["lp2"]

def th(mode, s0=0.0, dropq=True):
    cfg = dict(M.DEFAULT); cache = {}
    Sg = np.linspace(0.0005, 0.6, 60)
    vals = [lp2(S, cfg, mode, Rcache=cache, s0=s0, dropq=dropq) for S in Sg]
    k = int(np.argmin(vals)); a = Sg[max(k - 1, 0)]; c = Sg[min(k + 1, len(Sg) - 1)]
    r = minimize_scalar(lambda S: lp2(S, cfg, mode, s0=s0, dropq=dropq), bounds=(a, c), method="bounded",
                        options=dict(xatol=1e-7))
    Sb = r.x if r.fun < vals[k] else Sg[k]
    return lp2(Sb, cfg, mode, return_geom=True, s0=s0, dropq=dropq)

if __name__ == "__main__":
    for s0 in [0.0, 0.002]:
        sg, g = th("const", s0)
        print("const s0=%.3f: theta=%.6f" % (s0, sg), {k: round(v, 4) for k, v in g.items()})
    sg, g = th("zexp"); print("zexp (optimistic single-factor zeta expansion): theta=%.6f" % sg, {k: round(v, 4) for k, v in g.items()})
    sg, g = th("none"); print("none (Cauchy-Schwarz, P^2/Y dropped, slot-free): theta=%.6f" % sg, {k: round(v, 4) for k, v in g.items()})
