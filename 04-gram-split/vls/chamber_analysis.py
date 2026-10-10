"""Residue-adapted chamber analysis for the theta pairing T in three encodings:
   U  = (s, m, A)    exact after Patterson (A carries gamma_2), [[4,5,0],[5,0,1],[0,1,5]]      -- faithful, finite
   M0 = (t, A, m, s) A2_theory's 4-node encoding                                              -- conjugate cubic theta
   M1 = (t, A, m, s) corrected 4-node encoding [[1,2,0,0],[2,0,1,0],[0,1,0,5],[0,0,5,4]]       -- faithful, infinite type
For each chamber E (rows = new simple roots in old root coordinates) in which every residue root is +-simple:
  * free variables u'_k = <E_k, x> restricted to the residue plane (affine in (uA, um));
  * lengths l'_k = d(uA + M um)/du'_k   (Q = Z^M; units of log Z);
  * kappa_k = d um / du'_k : if all kappa_k <= 0 (not all 0) the chamber's convergence cone lies in Re um -> -infinity,
    so the contour must cross the coherent line P1 = {um = 1/3} whose residue is the Moebius signal S  ("P1 forced");
  * coefficient structure: residue numerators (Heckenberger m_kj), free-node Gauss types, free-free edge;
  * bound exponents beyond Q^{5/12} (random-size normalisation, sublattice factor K=1):
        trivial (l1+l2)/2 ; large sieve (bilinear only) ; CS target M/2.
"""
import sys
from math import gcd
sys.path.insert(0, '/Users/akashlevy/Documents/riemann/gram_split/vls')
from groupoid_bfs import bfs
from groupoid_search import Mof
import sympy as sp

uA, um, Ms = sp.symbols('uA um M')
MODELS = {
 'U':  dict(M=[[4, 5, 0], [5, 0, 1], [0, 1, 5]], names=['s', 'm', 'A'], res={'s': sp.Rational(1, 6)},
            free={'m': um, 'A': uA}),
 'M0': dict(M=[[5, 4, 0, 0], [4, 0, 1, 0], [0, 1, 0, 5], [0, 0, 5, 4]], names=['t', 'A', 'm', 's'],
            res={'t': sp.Rational(1, 3), 's': sp.Rational(1, 6)}, free={'A': uA, 'm': um}),
 'M1': dict(M=[[1, 2, 0, 0], [2, 0, 1, 0], [0, 1, 0, 5], [0, 0, 5, 4]], names=['t', 'A', 'm', 's'],
            res={'t': sp.Rational(1, 3), 's': sp.Rational(1, 6)}, free={'A': uA, 'm': um}),
}

def ls_exponent(l1, l2, e):
    """log_Z of the large-sieve constant for a bilinear form with kernel of order k = 6/gcd(6,e) (e != 0)."""
    k = 6 // gcd(6, e % 6)
    a, b = max(l1, 0.0), max(l2, 0.0)
    if k == 2:
        return max(a, b)
    if k == 3:   # Heath-Brown
        return max(a, b, 2 * (a + b) / 3)
    if k == 6:   # de Faveri Theta_6
        return max(a, b, 5 * a / 6 + b / 3, a / 3 + 5 * b / 6)
    return None

def analyse(name, cap=6000, maxcoef=30, maxlen=None, Mval=0.826, verbose=True):
    md = MODELS[name]
    dist, chi, mij = bfs(md['M'], cap=cap, maxcoef=maxcoef)
    names = md['names']; r = len(names)
    idx = {nm: i for i, nm in enumerate(names)}
    base = {nm: tuple(int(i == idx[nm]) for i in range(r)) for nm in names}
    reps = {}
    for E, w in dist.items():
        if maxlen is not None and len(w) > maxlen: continue
        kres = {}
        ok = True
        for nm in md['res']:
            ks = [k for k in range(r) if E[k] in (base[nm], tuple(-x for x in base[nm]))]
            if not ks: ok = False; break
            kres[nm] = ks[0]
        if not ok: continue
        free = [k for k in range(r) if k not in kres.values()]
        Mm = Mof(E, chi)
        x = sp.Matrix([md['res'].get(nm, md['free'].get(nm)) for nm in names])
        up = sp.Matrix(E) * x                     # new variables on the residue plane
        fv = [sp.expand(up[k]) for k in free]
        # invert: (uA, um) in terms of free new variables
        y1, y2 = sp.symbols('y1 y2')
        sol = sp.solve([sp.Eq(fv[0], y1), sp.Eq(fv[1], y2)], [uA, um], dict=True)
        if not sol:
            continue
        sol = sol[0]
        expo = sp.expand(sol[uA] + Ms * sol[um])
        lens = [sp.expand(sp.diff(expo, y)) for y in (y1, y2)]
        kap = [sp.diff(sol[um], y) for y in (y1, y2)]
        kapA = [sp.diff(sol[uA], y) for y in (y1, y2)]
        P1 = all(k <= 0 for k in kap) and any(k < 0 for k in kap)
        numer = {nm: {names_k: mij(E, kres[nm], kk) for kk, names_k in [(kk, kk) for kk in free]} for nm in kres}
        coupled = any(sum(1 for kk in free if mij(E, kres[nm], kk) > 0) >= 2 for nm in kres)
        resres = [Mm[kres[a]][kres[b]] for a in kres for b in kres if a < b]
        e = Mm[free[0]][free[1]]
        lnum = [float(l.subs(Ms, Mval)) for l in lens]
        triv = sum(max(l, 0) for l in lnum) / 2
        lsx = None if (coupled or e == 0) else ls_exponent(lnum[0], lnum[1], e) / 2
        key = (tuple(str(l) for l in lens), tuple(str(k) for k in kap), e, coupled,
               tuple((Mm[k][k]) for k in free), tuple(sorted((nm, tuple(sorted(v.items()))) for nm, v in numer.items())))
        prev = reps.get(key)
        if prev is None or len(w) < prev['wlen']:
            reps[key] = dict(w=''.join(names[i] for i in w) or 'id', wlen=len(w), lens=lens, kap=kap, kapA=kapA, P1=P1,
                             coupled=coupled, e=e, ftypes=[Mm[k][k] for k in free], numer=numer, resres=resres,
                             fv=fv, triv=triv, ls=lsx, lnum=lnum, mult=(prev['mult'] if prev else 0))
        reps[key]['mult'] += 1
    if verbose:
        print("=" * 110)
        print("model %s: %d chambers explored (cap %d, maxcoef %d), %d residue-adapted representation types"
              % (name, len(dist), cap, maxcoef, len(reps)))
        CS = Mval / 2
        for key, d in sorted(reps.items(), key=lambda kv: (kv[1]['P1'], max(kv[1]['lnum']))):
            sav = None
            best = min([v for v in (d['triv'], d['ls']) if v is not None])
            sav = CS - best
            print("-- %-22s mult=%-3d free lens=%s (at M=%.3f: %s)  um-coeffs kappa=%s  P1 forced=%s"
                  % (d['w'], d['mult'], [str(l) for l in d['lens']], Mval, [round(v, 3) for v in d['lnum']],
                     [str(k) for k in d['kap']], d['P1']))
            print("     free node types M_kk=%s edge=%d coupled=%s residue numerators=%s res-res edges=%s"
                  % (d['ftypes'], d['e'], d['coupled'], d['numer'], d['resres']))
            print("     bound exps (beyond Q^5/12): trivial=%.3f  LS=%s  CS=%.3f  => saving %.3f %s"
                  % (d['triv'], None if d['ls'] is None else round(d['ls'], 3), CS, sav,
                     '(but contains signal S)' if d['P1'] else ''))
    return reps

if __name__ == "__main__":
    which = sys.argv[1:] or ['U', 'M0', 'M1']
    for nm in which:
        analyse(nm, cap=(3000 if nm == 'M1' else 6000), maxcoef=(12 if nm == 'M1' else 30))
