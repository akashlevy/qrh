"""Redesign model (agent `redesign`): weighted high side with pluggable numerator bounds.

Units: everything inside a bin is in U-units (U = Z^d, d = h at the binding scale) unless noted.
E(bin) = base - l_y/2 + d*(R + sigma - z_c)   [weighted-type numerator: sum_C |S~_u| <= U^{R+sigma}]
E_old  = base - a*l_y + d*(R + delta/2 - z_c) [P78 pointwise numerator]
base   = a - B + h(z_c - 1/6) - ell/2 + q*ell, q = delta*x.
Counts R(delta,x) from gmodel.count (LP over witness splits, paper lemmas).
"""
import numpy as np
from gmodel import count

zc = 17/50
KAPPA = 0.75
CP = 1/(6*KAPPA)          # lem:plain prime-capacity rate (2/9)


def geom(ell, b, M=None):
    if M is None:
        M = 1-ell
    lx = (M-b)/2; ly = (M+b)/2; h = 1-lx+ell
    return lx, ly, h


def pieces(ell, b, B, dl, x, d=None, **kw):
    lx, ly, h = geom(ell, b)
    if d is None:
        d = h
    L = ell/d
    q = dl*x; a = (1+dl)/2
    Rv = count(dl, x, ell/h, **kw)
    base = a-B+h*(zc-1/6)-ell/2+q*ell
    n0 = 1-ly/d
    Eold = base-a*ly+d*(Rv+dl/2-zc)
    return dict(lx=lx, ly=ly, h=h, d=d, L=L, q=q, a=a, R=Rv, base=base, n0=n0, Eold=Eold)


def E_sigma(P, sigma):
    return P['base']-P['ly']/2+P['d']*(P['R']+sigma-zc)


def sigma_allowed(P):
    """largest sigma with E_sigma <= 0"""
    return -(P['base']-P['ly']/2)/P['d']-P['R']+zc


def sigma_wt(P, cP=CP):
    zN = min(P['L'], max(0.0, (1-2*P['n0'])*cP))
    return (1-P['R'])/4-P['q']*zN/2


# ---------------- designs ----------------
def des_paper(P):
    return min(P['Eold'], E_sigma(P, sigma_wt(P)))


def des_mollifier(theta_Z):
    """(a1) s-weight c = 1*lambda, lambda supported on q_d <= Z^theta_Z, |lambda_d| <= 1.
    Numerator = sum_d lambda_d chi_u(d) S_u(Y/q_d).  The only moment lemma (lem:plain) takes plain
    sums x prime polynomials, so the d-sum is removed by the triangle inequality:
    sum_d |lambda_d| q_d^{-1/2} ~ Z^{theta/2}; each plain sum has reflected length n0 + theta/d.
    The pointwise branch pays the same Z^{theta/2}."""
    def f(P):
        th = theta_Z/P['d']
        zN = min(P['L'], max(0.0, (1-2*(P['n0']+th))*CP))
        s = (1-P['R'])/4-P['q']*zN/2+th/2
        return min(P['Eold']+theta_Z/2, E_sigma(P, s))
    return f


def des_Lambda(P):
    """(a2) s-weight Lambda(s) (primes): numerator ~ -L'/L(w,chi_u); pointwise via the bin's
    zero-free region: |sum_{p~Y} chi_u(p)log p| << Y^{a+eps}; normalized: Y^{delta/2}."""
    s = P['ly']*(2*P['a']-1)/2/P['d']
    return min(P['Eold'], E_sigma(P, s))


def des_product(total_len='n0', cP=CP):
    """(b) numerator = S(n1) S(n2) (two plain sums, n1+n2 = t).  lem:plain bounds
    sum_u |S(n1)S(n2)Q|^2 << U with n1+n2+6kappa z <= 1: a SECOND moment of the numerator.
    Hoelder: sum_C |S1 S2| <= #C^{1/2} U^{1/2} U^{-q z}.
    Mixed Hoelder with |S_i|^4 moments is never better (dsigma/dalpha < 0, see NOTES)."""
    def f(P):
        t = P['n0'] if total_len == 'n0' else (P['ly']/P['d'] if total_len == 'Y' else total_len)
        z = min(P['L'], max(0.0, (1-t)*cP))
        s = (1-P['R'])/2-P['q']*z
        # pointwise branch for a product: Y^{-a} U^{delta} (two convexity factors)
        Eold2 = P['Eold']+P['d']*(2*P['a']-1)/2
        return min(Eold2, E_sigma(P, s))
    return f


def des_fractional_ideal(P):
    """Hypothetical: the 'ideal fractional moment' sum|S|^{2/n0} << U (no primes):
    sigma = (n0/2)(1-R).  This is the best a perfect product-splitting could approach if the
    needed higher diagonal moments existed (they do not: sextic LS diagonal only to length 1/2)."""
    s = P['n0']/2*(1-P['R'])
    return min(P['Eold'], E_sigma(P, s), E_sigma(P, sigma_wt(P)))


def des_sparseLS(P):
    """Hypothesis (sLS): sum_{u in C}|S~_u|^2 << U^eps (#C + Y/1) i.e. Bessel for the bad-row vectors.
    sigma = max(0,(ly/d - R)/2)."""
    s = max(0.0, (P['ly']/P['d']-P['R'])/2)
    return min(P['Eold'], E_sigma(P, s), E_sigma(P, sigma_wt(P)))


def des_annihilate(P, sstar_floor=None):
    """(c) optimistic target-dependent annihilation: if the number of linear conditions
    #C * D fits in the Y-dimensional s-space (R + theta_D < ly/d), the bin is declared closed;
    theta_D = truncation of the Euler-correction twists k <= U^theta_D needed to bring the
    untwisted remainder down to sigma*: theta_D = (sigma_wt - sigma*)/(delta/2).
    Costs of c1 on other bins and on the low side are IGNORED (optimistic)."""
    sw = sigma_wt(P); ss = sigma_allowed(P)
    if min(P['Eold'], E_sigma(P, sw)) <= 0:
        return min(P['Eold'], E_sigma(P, sw))
    dl = 2*P['a']-1
    thD = max(0.0, (sw-ss)/(dl/2))
    if P['R']+thD < P['ly']/P['d']:
        return -1e-9
    return min(P['Eold'], E_sigma(P, sw))


def worst(design, ell, b, B, grid_d=np.linspace(0.05, 0.75, 15), grid_x=(0, 0.25, 0.5), **kw):
    w = -9; arg = None
    for dl in grid_d:
        for x in grid_x:
            P = pieces(ell, b, B, dl, x, **kw)
            e = design(P)
            if e > w:
                w = e; arg = (round(dl, 3), x, round(P['R'], 4))
    return w, arg


def best_b(design, ell=0.1773, B=41/47, bs=np.linspace(0.0, 0.2, 11), **kw):
    return min((worst(design, ell, b, B, **kw)+(b,)) for b in bs)


if __name__ == "__main__":
    ell = 0.1773; B = 41/47
    designs = [('paper (weighted)', des_paper),
               ('(a1) mollifier theta=0.005', des_mollifier(0.005)),
               ('(a1) mollifier theta=0.02', des_mollifier(0.02)),
               ('(a2) Lambda s-weight', des_Lambda),
               ('(b) product t=n0', des_product('n0')),
               ('(b) product t=n0, cP=1', des_product('n0', cP=1.0)),
               ('(b) product Y=Y1Y2 unrefl', des_product('Y')),
               ('HYP ideal fractional mom.', des_fractional_ideal),
               ('(c) annihilation (optim.)', des_annihilate),
               ('HYP sparse LS (Bessel on C)', des_sparseLS)]
    for name, f in designs:
        w, arg, b = best_b(f, ell, B)
        print(f"{name:30s} worst={w:+.5f} at (delta,x,R)={arg} b={b:.3f}")
