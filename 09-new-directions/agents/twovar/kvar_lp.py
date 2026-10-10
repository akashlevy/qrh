"""
k-Poisson-variable probes in P78 Part II shape (agent `twovar`, 10 Oct 2026).  MODEL.

Probe (k = number of Poisson-transformed variables, Jacobi triple split as k Poisson + row):

    I = prod_i Q_i^{-1/2}  sum_{m_1..m_k}  prod_i A^{(i)}_{m_i}  B_{m_1 ... m_k},

    A^{(i)}_m = Y_i^{-3/2} sum_{s ~ Y_i} c_i(s) g(psi_s^{e_i}, -m)      (Gauss-sum polynomial, Q_i = X_i Y_i)
    B_n       = sum_{A ~ Z} b(A) q_A^{-1/2} kappa_A(n)  [x prime marks 1_{p|A}, Z -> Zp]

  k = 1 cubic : b = gamma_1^2 (theta(x)theta, c_B = 4 by cB-theory; alt-detect template)
  k = 1 sextic: b = cubic theta, sextic twist -> quadratic residual (P78)
  k = 2 cubic : b = single cubic theta (gamma_1), twist kappa_A(m1 m2), full index, c_B = 2
  k = 3 cubic : b = lambda-bar eta (Hecke), degree 1, c_B = 1

Exponents in units of log Z: X_i = Z^{l_i}, Y_i = Z^{y_i}, P = Z^ell, M_i = l_i + y_i,
F_i = Z P / X_i (Poisson frequency range of variable i, f_i = 1 + ell - l_i).

Derived pieces (see NOTES.md sec 2 for the derivation):
  signal   S = prod_i F_i^{c_i - 1/2} * P^{1/2} * Z^{beta - 1/2}      (copies prod F_i^{c_i})
  low      L = prod_i X_i^{1/2} G_i^{1/2} * E^{1/2}
           G_i = 1 + P_{a,i}^{gb_i} + gq P_{a,i}^2/Y_i ,  P_{a,i} = Y_i/X_i     (Gram of A^{(i)})
           E   = energy excess of sum_{m ~ prod Q_i} d_k(m)|B_m|^2 (family M' = sum M_i, row length 1+l')
  rows     (GRH for rows) for every subset J of variables whose frequency is a copy (u_i = cube):
           R_J = prod_{i in J} F_i^{c_i - 1/2} * prod_{i not in J} (Z P / Q_i)^{1/2}      (J != all)
  compensation (d = ell): unmarked probe at (X_i/P, Y_i/P), relative weight P^{-comp_drop}.

For k = 1 these are exactly new-architectures/template_lp.py and alt-detect/p78_cubic_lp.py
(validated in validate.py: sextic 13/15, cubic c_B = 2/3/4 -> 13/17, 0.8596, 10/11).
"""
import itertools, json, sys
import numpy as np
from scipy.optimize import linprog


def default_cfg(kind="cubic2"):
    base = dict(k=2, c=[1/3, 1/3], gb=None, gq=0.0, cB=2, ls=3, caseA=None, third=True, wE=0.0,
                margin=0.05, margin_y=None, comp=True, comp_drop=None, marks=True, partial=True,
                lmax=2.0, ymax=6.0, ellmax=2.0, Fge1=False, fix=None, sym=False, extra=None)
    presets = {
        "sextic1": dict(k=1, c=[1/6], gq=1.0, cB=2, ls=2, caseA="paper", third=True),
        "cubic1": dict(k=1, c=[1/3], gq=0.0, cB=2, ls=3, caseA=None, third=True),
        "cubic2": dict(k=2, c=[1/3, 1/3], cB=2, ls=3, third=True, Fge1=True),
        "cubic3": dict(k=3, c=[1/3, 1/3, 1/3], cB=1, ls=3, third=False, Fge1=True),
    }
    base.update(presets[kind])
    return base


def solve(cfg, verbose=False):
    k = cfg["k"]; cs = list(cfg["c"])
    gbs = list(cfg["gb"]) if cfg.get("gb") is not None else list(cs)
    gq, cB, ls = cfg["gq"], cfg["cB"], cfg["ls"]
    nv = 2 * k + 1                      # l_1..l_k, y_1..y_k, ell
    def e(j):
        v = np.zeros(nv + 1); v[j] = 1.0; return v
    ONE = np.zeros(nv + 1); ONE[-1] = 1.0
    L = [e(i) for i in range(k)]; Yv = [e(k + i) for i in range(k)]; ELL = e(2 * k)
    Z0 = np.zeros(nv + 1)
    Mi = [L[i] + Yv[i] for i in range(k)]
    Mtot = sum(Mi)
    comp_drop = cfg["comp_drop"] if cfg.get("comp_drop") is not None else (k + 1) / 2.0
    wE = cfg["wE"]

    def nodes(Mp, lp):
        f = Mp; l = ONE + lp
        lR = cB * (Mp + lp) - (ONE + lp)
        if ls == 2:
            nd = [[Z0, lR - f - wE * ONE]]
        else:
            nd = [[Z0, l - f, (2 / 3) * (l + f) - f],
                  [Z0, l - (2 / 3) * f],
                  [Z0, lR - f - wE * ONE, (2 / 3) * (lR + f) - f - wE * ONE],
                  [Z0, lR - (2 / 3) * f - wE * ONE]]
        if cfg.get("Eperfect"):
            nd = [[Z0]]
        if cfg.get("node_filter") is not None:
            nd = [nd[j] for j in cfg["node_filter"]]
        return nd

    lams = (0.0, 1.0) if cfg["comp"] else (0.0,)
    nn = len(nodes(Mtot, ELL))
    best = None
    for choice in itertools.product(range(nn), repeat=len(lams)):
        A, B = [], []
        def ge(form):                     # sigma >= form
            A.append(list(form[:nv]) + [-1.0]); B.append(-form[nv])
        def le0(form):                    # form <= 0
            A.append(list(form[:nv]) + [0.0]); B.append(-form[nv])
        for lam, ch in zip(lams, choice):
            d = lam * ELL
            Mp = Mtot - lam * 2 * k * ELL
            lp = (1 - lam) * ELL
            base = 0.5 * ONE - 0.5 * ELL - lam * comp_drop * ELL
            for i in range(k):
                base = base + (0.5 - cs[i]) * (ONE + ELL) + cs[i] * L[i]
            gp_i = []
            for i in range(k):
                b = Yv[i] - L[i]
                g = [Z0]
                if gbs[i] > 0: g.append(gbs[i] * b)
                if gq > 0: g.append(gq * (2 * b - Yv[i] + d))
                gp_i.append(g)
            ep = list(nodes(Mp, lp)[ch])
            if cfg["caseA"] == "paper":
                ep.append(0.75 * lp - 0.5 * Mp + 0.25 * ONE)
            if cfg["third"]:
                ep.append((1 / 3) * Mp + (4 / 3) * lp - (2 / 3) * ONE - wE * ONE)
            for gs in itertools.product(*gp_i):
                for en in ep:
                    ge(base + 0.5 * sum(gs) + 0.5 * en)
        # ---- high side, GRH for rows: every proper subset J of 'copy' variables
        subsets = [()] if not cfg["partial"] else [J for r in range(k) for J in itertools.combinations(range(k), r)]
        for J in subsets:
            form = 0.5 * ONE - 0.5 * ELL
            for i in range(k):
                if i in J: continue
                form = form + 0.5 * (ONE + ELL - Mi[i]) + (0.5 - cs[i]) * (ONE + ELL - L[i])
            ge(form)
        # ---- side conditions (template: M - 2 ell >= m, ly - ell >= m, lx >= ell)
        m = cfg["margin"]; my = cfg["margin_y"] if cfg.get("margin_y") is not None else m
        if cfg["comp"]:
            le0(-(Mtot - 2 * k * ELL) + m * ONE)
        for i in range(k):
            le0(-Yv[i] + ELL + my * ONE)
            le0(-L[i] + ELL)
            if cfg["Fge1"]:
                le0(L[i] - ONE - ELL)      # F_i >= 1
        if not cfg["marks"]:
            le0(ELL)
        if cfg.get("sym") and k >= 2:
            for i in range(1, k):
                for s in (1, -1):
                    le0(s * (L[i] - L[0])); le0(s * (Yv[i] - Yv[0]))
        if cfg.get("fix"):
            for name, val in cfg["fix"].items():   # e.g. {"y2": 0.0}
                kind, idx = name[0], (int(name[1:]) - 1 if len(name) > 1 and name[1:].isdigit() else None)
                var = {"l": L, "y": Yv}.get(kind)
                f = (var[idx] if var is not None else ELL) - val * ONE
                le0(f); le0(-f)
        if cfg.get("extra"):
            for fn in cfg["extra"]:
                fn(ge, le0, L, Yv, ELL, ONE)
        bounds = [(0, cfg["lmax"])] * k + [(0, cfg["ymax"])] * k + [(0, cfg["ellmax"]), (None, None)]
        cost = np.zeros(nv + 1); cost[-1] = 1.0
        res = linprog(cost, A_ub=np.array(A), b_ub=np.array(B), bounds=bounds, method="highs")
        if res.status != 0:
            continue
        if best is None or res.x[-1] < best["theta"] - 1e-12:
            x = res.x
            best = dict(theta=float(x[-1]), l=[float(v) for v in x[:k]], y=[float(v) for v in x[k:2 * k]],
                        ell=float(x[2 * k]), M=float(sum(x[:2 * k])), nodes=choice)
    return best


def binding(cfg, sol, tol=1e-7):
    """recompute all constraints at the optimum and list the tight ones (for reporting)."""
    k = cfg["k"]; cs = cfg["c"]; gbs = cfg["gb"] if cfg.get("gb") is not None else cs
    l = np.array(sol["l"]); y = np.array(sol["y"]); ell = sol["ell"]; th = sol["theta"]
    out = []
    # low (d=0) pieces
    base = 0.5 - ell / 2 + sum((0.5 - cs[i]) * (1 + ell) + cs[i] * l[i] for i in range(k))
    g = sum(max(0, gbs[i] * (y[i] - l[i]), cfg["gq"] * (2 * (y[i] - l[i]) - y[i])) for i in range(k))
    out.append(("low base (signal/low, no Gram/energy)", base))
    out.append(("low + Gram/2", base + g / 2))
    Mt = l.sum() + y.sum()
    for J in [J for r in range(k) for J in itertools.combinations(range(k), r)]:
        v = 0.5 - ell / 2 + sum(0.5 * (1 + ell - l[i] - y[i]) + (0.5 - cs[i]) * (1 + ell - l[i]) for i in range(k) if i not in J)
        out.append((f"GRH rows, copy set J={J}", v))
    return [(n, v, abs(v - th) < 1e-6) for n, v in out]


if __name__ == "__main__":
    kind = sys.argv[1] if len(sys.argv) > 1 else "cubic2"
    over = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    cfg = default_cfg(kind); cfg.update(over)
    print(json.dumps(solve(cfg)))
