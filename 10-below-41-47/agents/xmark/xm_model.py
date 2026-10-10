"""xmark: payoff of a cross-mark saving in lem:marked (weighted-numerator high side, count LP).

Marked capacity models (U-units, m = 1):
  * absolute saving delta (energy exponent saved on the cross-tuple terms):
        excess(r,z) = max(r+2z-1, cubic(r,z)) <= delta,
    cubic branch options: 'keep'  -> 2r+8z <= 3          (cap not improved)
                          'relax' -> (2r+8z-3)/4 <= delta (cap improved by the same energy saving)
                          'drop'  -> no cubic branch
  * relative saving phi (fraction of the converted mark length): z <= (1-r)/(2-phi)  [mk = 1/(2-phi)]
Supply: physical L = ell/h, or unlimited (polarised slot powers, amplitude agent Prop 3.3).
Plain witness: capacity (1-2m)*cP (paper cP = 2/9); long witness slope alpha (paper 5/6).
High side: E_old, E_wt exactly as in model/gmodel.py (kappa = 3/4, z_c = 17/50)."""
import numpy as np
zc = 17/50

def zm_fun(r, delta=0.0, cubic='keep', phi=None):
    if phi is not None:
        z = (1-r)/(2-phi)
    else:
        z = (1+delta-r)/2
    if cubic == 'keep':
        z = np.minimum(z, (3-2*r)/8)
    elif cubic == 'relax':
        z = np.minimum(z, (3+4*delta-2*r)/8)
    return np.maximum(z, 0)

def count(dl, x, L, delta=0.0, cubic='keep', phi=None, cP=2/9, alpha=5/6, nt=101, nr=201):
    q = dl*x
    T = np.linspace(1, 1.5, nt)[:, None]; s = np.linspace(0, 1, nr)[None, :]
    r = np.maximum(0, T-0.5) + (T-np.maximum(0, T-0.5))*s; m = T-r
    zp = np.clip((1-2*m)*cP, 0, L); b = 1-dl*2*m-2*q*zp
    zm = np.minimum(L, zm_fun(r, delta, cubic, phi))
    a = np.where(r < 1, 1-dl*r-2*q*zm, 1-dl+(alpha-dl)*(r-1))
    v = np.minimum(a, b)
    return max(float(np.min(np.max(v, axis=1))), 1-dl)

def E(ell, b, B, dl, x, kappa=0.75, unlimited=False, **kw):
    M = 1-ell; lx = (M-b)/2; ly = (M+b)/2; h = 1-lx+ell; d = h; L = ell/h
    q = dl*x; a = (1+dl)/2
    Rv = count(dl, x, 10.0 if unlimited else L, **kw)
    base = a-B+h*(zc-1/6)-ell/2+q*ell
    Eold = base-a*ly+d*(Rv+dl/2-zc)
    n0 = 1-ly/d; zN = min(L, max(0, (1-2*n0)/(6*kappa)))
    Ewt = base-ly/2+d*(Rv+(1-Rv)/4-q*zN/2-zc)
    return min(Eold, Ewt), Rv

def worst(ell=0.1773, B=41/47, bs=np.linspace(0, 0.10, 11), gd=np.linspace(0.2, 0.75, 56),
          gx=(0, 0.2, 0.3, 0.4, 0.45, 0.5), **kw):
    best = None
    for b in bs:
        w = -9; arg = None
        for dl in gd:
            for x in gx:
                e, Rv = E(ell, b, B, dl, x, **kw)
                if e > w: w = e; arg = (round(float(dl), 3), x, round(Rv, 4))
        if best is None or w < best[0]: best = (w, arg, round(float(b), 3))
    return best

if __name__ == "__main__":
    import sys, json
    rows = []
    rows.append(('paper (delta=0, cap kept, physical supply)', dict()))
    for dlt in (0.02, 0.05, 0.10, 0.13, 0.15, 0.2):
        for cub in ('keep', 'relax', 'drop'):
            for unl in (False, True):
                rows.append((f'delta={dlt} cubic={cub} unlimited={unl}', dict(delta=dlt, cubic=cub, unlimited=unl)))
    for phi in (0.3, 0.5, 0.6, 0.656, 0.7, 1.0):
        rows.append((f'phi={phi} (mk={1/(2-phi):.3f}) cubic=drop unlimited', dict(phi=phi, cubic='drop', unlimited=True)))
    for lab, kw in rows:
        w, arg, b = worst(**kw)
        print(json.dumps(dict(cfg=lab, worst=round(w, 5), delta_x_R=arg, b=b)), flush=True)
