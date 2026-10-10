# Sanity check: the reflection row scalar for a prime row R = p contains
#   alpha(p)^{-2} * tau^-_{p,1} * tau^+_{p,3}   (times symbols of fixed elements, i.e. ray-class phases)
# Is this a ray-class function of p?  Compute w(p) and test constancy within classes mod 72.
import sys, os, numpy as np
os.chdir("/Users/akashlevy/Documents/riemann/gram_split"); sys.path.insert(0, ".")
from zomega import *
def gam(Pr, j, sign=+1):
    p = Pr.q; x = np.arange(1, p, dtype=np.int64); jj = Pr.symbol(x, np.zeros_like(x))
    ph = np.exp(sign * -2j * np.pi * ((Pr.b * x) % p) / p)
    return np.sum(Z6[(j * jj) % 6] * ph) / np.sqrt(Pr.N)
P = [Pr for Pr in primes_upto(60000) if Pr.kind == "split"]
rows = []
for Pr in P:
    al = cval(Pr.a, Pr.b); al = al / abs(al)
    t1m = gam(Pr, 1, sign=-1)    # tau^-_{p,1}
    t3p = gam(Pr, 3, sign=+1)    # tau^+_{p,3}
    g2 = gam(Pr, 2, sign=+1)
    w = np.conj(al) ** 2 * t1m * t3p
    rows.append(((Pr.a % 72, Pr.b % 72), w, al, g2, t3p))
from collections import defaultdict
cls = defaultdict(list)
for k, w, al, g2, t3 in rows: cls[k].append((w, al, g2, t3))
big = sorted(cls.items(), key=lambda kv: -len(kv[1]))[:5]
for k, L in big:
    W = np.array([x[0] for x in L]); ang = np.angle(W)
    # how many distinct values (rounded) and circular spread
    distinct = len(set(np.round(ang, 3)))
    R = abs(np.mean(W)); R6 = abs(np.mean(W ** 6)); R12 = abs(np.mean(W**12))
    # DH prediction: w = (ray) * conj(gamma_2)^4  ->  w * gamma_2^4 should take finitely many values
    v = np.array([x[0] * x[2] ** 4 for x in L])
    dv = len(set(np.round(np.angle(v), 3)))
    print("class mod 72 %s: #primes=%d  distinct arg(w)=%d  |mean w|=%.3f |mean w^6|=%.3f |mean w^12|=%.3f   distinct arg(w*gamma2^4)=%d"
          % (k, len(L), distinct, R, R6, R12, dv))
# all classes: is w*gamma2^4 constant per class mod 72, while w itself is not?
nc = 0; bad = 0; spread = []
for k, L in cls.items():
    if len(L) < 2: continue
    nc += 1
    v = np.array([x[0] * x[2] ** 4 for x in L]); W = np.array([x[0] for x in L])
    if np.max(np.abs(v - v[0])) > 1e-6: bad += 1
    spread.append(np.max(np.abs(W - W[0])))
print("classes with >=2 primes: %d ; classes where w*gamma2^4 not constant: %d ; median max|w-w0| within class: %.3f" % (nc, bad, float(np.median(spread))))
print("total split primes:", len(rows))
