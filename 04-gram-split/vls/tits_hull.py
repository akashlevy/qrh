"""C.1 numerics: the reflection group of the second-moment series
   Z2(z1,z2,u) = sum_mu N(mu)^{-u} R(z1;mu) conj R(z2;mu),
in centred unitary coordinates a = Re z1' - 1/2, b = Re z2' - 1/2, c = Re u - 1/2.
Reflections:
  s_c : (a,b,c) -> (a+c, b+c, -c)        mu-Poisson (Hecke L(u,chi_{m1 m2})), conductor N(m1 m2)   [R]
  s_a : (a,b,c) -> (-a, b, c + k a)      functional equation of R(z1;mu), conductor N(mu)^k       [H], k = 3 (A2_theory)
  s_b : (a,b,c) -> (a, -b, c + k b)
Absolute-convergence region C0 = {a,b,c > 1/2}.  The bound obtained at a point is
  sum_mu |F'(mu)|^2 << Q^{5/6} Q^{a+b} L^{1/2+c},   L = Q^lam,
so the exponent (in Q-units, after removing Q^{5/6}) is f = a + b + lam (c + 1/2).
Large sieve = 1 (L<Q); target (V'_LS) = f < 1.
We (1) print the Cartan matrix / Coxeter data, (2) enumerate the group to length N,
(3) minimise f over apexes g(1/2,1/2,1/2) of the image cones g(C0) and report whether any image
cone has a recession direction with f_lin < 0, (4) the same restricted to images that do not cross
the polar hyperplanes c = 1/2 (u=1, diagonal) and its images (crude: report the f-value of the
apex of every image cone together with whether the segment from C0's apex to it crosses a polar
hyperplane of the list)."""
import numpy as np, itertools, sys
lam = float(sys.argv[1]) if len(sys.argv) > 1 else 0.79
k = int(sys.argv[2]) if len(sys.argv) > 2 else 3
Sa = np.array([[-1,0,0],[0,1,0],[k,0,1]],float)
Sb = np.array([[1,0,0],[0,-1,0],[0,k,1]],float)
Sc = np.array([[1,0,1],[0,1,1],[0,0,-1]],float)
gens = {'a':Sa,'b':Sb,'c':Sc}
# Cartan / Coxeter check: products of pairs, orders
def order(M, cap=50):
    P = np.eye(3)
    for n in range(1, cap+1):
        P = P @ M
        if np.allclose(P, np.eye(3)): return n
    return None
for x,y in [('a','c'),('b','c'),('a','b')]:
    print("order of s_%s s_%s :"%(x,y), order(gens[x]@gens[y]))
# invariant bilinear form?  B with S^T B S = B for all generators
import scipy.linalg as sl
A = []
for S in gens.values():
    # S^T B S - B = 0, B symmetric: 6 unknowns
    rows = []
    idx = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    for i in range(3):
        for j in range(i,3):
            row = []
            for (p,q) in idx:
                E = np.zeros((3,3)); E[p,q]=1; E[q,p]=1
                D = S.T@E@S - E
                row.append(D[i,j])
            rows.append(row)
    A += rows
A = np.array(A); ns = sl.null_space(A)
print("dim of invariant symmetric forms:", ns.shape[1])
if ns.shape[1]:
    v = ns[:,0]; idx = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    B = np.zeros((3,3))
    for t,(p,q) in enumerate(idx):
        B[p,q] += v[t]; 
        if p!=q: B[q,p] += v[t]
    B = B/np.max(np.abs(B))
    print("invariant form B =\n", np.round(B,4)); print("signature eigenvalues:", np.round(np.linalg.eigvalsh(B),4))
# enumerate group
N = int(sys.argv[3]) if len(sys.argv) > 3 else 14
elems = {(): np.eye(3)}
frontier = [()]
for L in range(N):
    nf = []
    for w in frontier:
        for g,S in gens.items():
            if w and w[-1]==g: continue
            M = S @ elems[w]
            key = w+(g,)
            # dedupe by matrix
            elems[key] = M; nf.append(key)
    frontier = nf
mats = {}
for w,M in elems.items():
    kk = tuple(np.round(M,6).ravel())
    if kk not in mats or len(w) < len(mats[kk][0]): mats[kk] = (w,M)
print("distinct group elements up to length %d: %d" % (N, len(mats)))
flin = np.array([1.0,1.0,lam])
apex0 = np.array([0.5,0.5,0.5])
best = []; neg_dirs = 0
for kk,(w,M) in mats.items():
    ap = M@apex0
    val = flin@ap + lam/2
    dirs = [flin@(M@e) for e in np.eye(3)]
    if min(dirs) < -1e-12: neg_dirs += 1
    best.append((val, w, ap, min(dirs)))
best.sort(key=lambda t: t[0])
print("lam=%.3f k=%d ; # image cones with a recession direction of negative f: %d" % (lam,k,neg_dirs))
print("smallest apex values f (need <1 for (V'_LS); large sieve =1):")
for val,w,ap,md in best[:12]:
    print("  f=%.4f  word=%s  apex=%s  min f_lin on generators=%.3f" % (val, ''.join(w) or 'id', np.round(ap,3), md))
