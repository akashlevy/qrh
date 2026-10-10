"""Checks of cubic_data.c output against brute force (small norms).
 (1) G~(pi) vs direct sum over O/pi (split and inert), |G~| = 1, G~^3 = -pi/|pi|
 (2) twisted multiplicativity for composite squarefree n vs direct sum over O/n (O/n enumerated as
     u + v w, 0 <= u,v < N(n), divided by N(n) for multiplicity)
 (3) cubic reciprocity (h/n)_3 = (n/h)_3 for primary h, n (symbols via prime factorisations)
usage: python3 verify.py <prefix>"""
import sys, itertools, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__file__))
import cbn_lib as L

pref = sys.argv[1]
P, S = L.load(pref)
print(f"{len(P)} primes, {len(S)} squarefree n, max N = {S['N'].max()}")
mods = {}
def mod(i):
    if i not in mods: mods[i] = L.Mod(P[i])
    return mods[i]

def sym_n(A, B, idxs):
    """(x / n)_3 for x = A + B w (arrays), n = prod P[idxs]; exponent or -1"""
    return L.chi_exp([mod(i) for i in idxs], np.atleast_1d(A), np.atleast_1d(B))

def direct_gauss(a, b, idxs):
    Nn = int(L.norm(a, b))
    u, v = np.meshgrid(np.arange(Nn), np.arange(Nn), indexing="ij"); u = u.ravel(); v = v.ravel()
    j = sym_n(u, v, idxs)
    ok = j >= 0
    # Tr(x/n) = Tr(x conj(n))/N(n); x conj(n) = (u + v w)(a' + b' w), conj(n) = (a-b) - b w
    ca, cb = a - b, -b
    xa, xb = L.mulw(u, v, ca, cb)
    tr = (2 * xa - xb) % Nn          # Tr(xa + xb w) = 2 xa - xb
    s = np.sum(L.W3 ** j[ok] * np.exp(2j * np.pi * tr[ok] / Nn)) / Nn   # multiplicity Nn
    return s / np.sqrt(Nn)

# (1) primes
G = P["gre"] + 1j * P["gim"]
print("max ||G|-1| =", np.max(np.abs(np.abs(G) - 1)))
pi = L.cval(P["a"], P["b"])
print("max |G^3 + pi/|pi|| =", np.max(np.abs(G ** 3 + pi / np.abs(pi))))
errs = []
for i in range(min(40, len(P))):
    if P["N"][i] > 1500: break
    d = direct_gauss(int(P["a"][i]), int(P["b"][i]), [i])
    errs.append(abs(d - G[i]))
print(f"direct vs C, first {len(errs)} primes (N<=1500): max err {max(errs):.2e}")

# (2) composites: find squarefree n among S with 2-3 prime factors and N <= 900
pidx = {(int(P['a'][i]), int(P['b'][i])): i for i in range(len(P))}
small = [i for i in range(len(P)) if P["N"][i] <= 100]
cnt = 0; mx = 0
for k in (2, 3):
    for comb in itertools.combinations(small, k):
        a, b = 1, 0; Nn = 1
        for i in comb:
            a, b = L.mulw(a, b, int(P['a'][i]), int(P['b'][i])); Nn *= int(P['N'][i])
        if Nn > 1000: continue
        sel = np.nonzero((S['A'] == a) & (S['B'] == b))[0]
        assert len(sel) == 1, (a, b, comb)
        s = S[sel[0]]
        Gc = s['gre'] + 1j * s['gim']
        d = direct_gauss(a, b, list(comb))
        assert s['mu'] == (-1) ** k
        mx = max(mx, abs(d - Gc)); cnt += 1
        if cnt <= 6:
            print(f"  n={a:+d}{b:+d}w N={Nn} factors N={[int(P['N'][i]) for i in comb]}: C {Gc:.6f} direct {d:.6f}")
print(f"twisted multiplicativity: {cnt} composites, max err {mx:.2e}")

# (3) reciprocity
rng = np.random.default_rng(1)
bad = 0; tot = 0
idx_small = [i for i in range(len(P)) if P['N'][i] < 3000]
for _ in range(300):
    i1, i2 = rng.choice(idx_small, 2, replace=False)
    if P['p'][i1] == P['p'][i2] and P['kind'][i1] == P['kind'][i2]:
        continue
    j12 = mod(i1)(np.array([P['a'][i2]]), np.array([P['b'][i2]]))[0]   # (pi2 / pi1)
    j21 = mod(i2)(np.array([P['a'][i1]]), np.array([P['b'][i1]]))[0]   # (pi1 / pi2)
    tot += 1; bad += int(j12 != j21)
print(f"cubic reciprocity (pi1/pi2) = (pi2/pi1): {tot - bad}/{tot} agree")
