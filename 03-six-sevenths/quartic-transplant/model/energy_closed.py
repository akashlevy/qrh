"""Closed forms of the energy loss e(M', l') (log_Z exponent above the row count M'), verified against energy.py LPs."""
import numpy as np

def e_cubic_paper(Mp, lp, w=0.0):
    """7/8 paper (Lemma probe-row-norm / reflected-energy), paper count O/2; case A only when a dual dyad exists."""
    return max(0.0, ((1 + 3 * lp - 2 * Mp) / 4) if 2 * Mp + lp >= 1 else 0.0, Mp + lp - 1 - w, (Mp + 4 * lp - 2) / 3 - w)

def _q_reflected(Mp, lp, pieces):
    T0 = 3 * Mp + 2 * lp - 1                       # dual length: rows h^3 (3H) + marks p^2 & length (3 z_a) - N_*, N_* = 1+l'
    if T0 < 0:
        return 0.0                                  # conductor below length: completed sum negligible (no polar term)
    vals = [0.0, lp - T0 / 4]                       # short-dual branch: uncancelled mark sum, powerful-part trade 1/4
    for sH, sv in pieces:                           # full dual v = T0, H = M', z_a = l'
        vals.append(sH * Mp + sv * T0 - Mp)
    return max(vals)

DF4 = [(1, 0), (0, 1), (0.75, 0.5), (0.5, 0.75)]
PERF = [(1, 0), (0, 1)]
LOA = [(1, 0)]
BGL4 = [(1, 0), (0, 1), (2 / 3, 2 / 3)]

def e_quartic_reflected(Mp, lp, sieve="dF"):
    return _q_reflected(Mp, lp, {"dF": DF4, "perfect": PERF, "LoA": LOA, "BGL": BGL4}[sieve])

def e_quartic_direct(Mp, lp, sieve="dF"):
    """unreflected: rows M', columns N_* = 1 + l' (completed Kubota sum incl. marked slots)."""
    N = 1 + lp
    pcs = {"dF": DF4, "perfect": PERF, "LoA": LOA, "BGL": BGL4}[sieve]
    return max(sH * Mp + sv * N for sH, sv in pcs) - Mp

def e_quartic(Mp, lp, sieve="dF", allow_direct=True):
    r = e_quartic_reflected(Mp, lp, sieve)
    if allow_direct:
        r = min(r, max(0.0, e_quartic_direct(Mp, lp, sieve)))
    return r

if __name__ == "__main__":
    import sys; sys.path.insert(0, ".")
    from energy import energy_loss, energy_loss_refined, CUBIC, QUARTIC_DF, QUARTIC_PERFECT, QUARTIC_LOA, QUARTIC_BGL
    rng = np.random.default_rng(11)
    tests = [("dF", QUARTIC_DF), ("perfect", QUARTIC_PERFECT), ("LoA", QUARTIC_LOA), ("BGL", QUARTIC_BGL)]
    for name, cfg in tests:
        worst = 0; arg = None
        for _ in range(2500):
            Mp = rng.uniform(0.05, 1.3); lp = rng.uniform(0, 0.5)
            a = energy_loss_refined(Mp, lp, cfg); b = e_quartic_reflected(Mp, lp, name)
            if abs(a - b) > worst: worst = abs(a - b); arg = (round(Mp, 4), round(lp, 4), a, b)
        print(f"quartic {name:8s} refined-count LP vs closed form: max|diff| = {worst:.2e}  {arg}")
    worst = 0
    for _ in range(2500):
        Mp = rng.uniform(0.05, 1.3); lp = rng.uniform(0, 0.5)
        a = energy_loss(Mp, lp, CUBIC); b = e_cubic_paper(Mp, lp)
        T = 2 * Mp + lp - 1
        if T >= 0: worst = max(worst, abs(a - b))
    print(f"cubic paper-count LP vs closed form (dual nonempty): max|diff| = {worst:.2e}")
