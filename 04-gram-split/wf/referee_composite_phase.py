# Referee check: for squarefree m = p1 p2 (split primes over distinct rational primes), build the
# m-dependent part of the reflection row scalar exactly as in paper78 eq. (reflection-row-phase), c_F = 1:
#   W(m) = conj(alpha(m))^2 * prod_{p|m} chi_p(sigma_p)^{-2} tau^-_{p,1} tau^+_{p,3} chi_p(eps_p)^{-3},
#   sigma_p = lambda^2 m/p,  eps_p = -((lambda^3 m/p) sigma_p)^{-1}.
# Test (a) W(m) * gamma_2(m)^4 is constant on classes m mod 72 (=> ray class), with gamma_2(m) the
# composite Gauss sum computed directly mod m; (b) W(m) itself is not.
import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
from collections import defaultdict
P = [Pr for Pr in primes_upto(700) if Pr.kind == "split" and Pr.N > 7]
def gam_p(Pr, j, sign):
    p = Pr.q; x = np.arange(1, p, dtype=np.int64); jj = Pr.symbol(x, np.zeros_like(x))
    return np.sum(Z6[(j * jj) % 6] * np.exp(sign * -2j * np.pi * ((Pr.b * x) % p) / p)) / np.sqrt(Pr.N)
def sym(Pr, a, b):
    return int(Pr.symbol(np.array([a]), np.array([b]))[0])
def inv_mod_p(Pr, a, b):
    # inverse of a + b w modulo split prime Pr, as an integer residue mod p (w = r mod p)
    p = Pr.q; v = (a + b * Pr.r) % p; return pow(int(v), -1, p)
lam = (1, 2)   # lambda = 1 + 2w
def powel(a, b, k):
    x = (1, 0)
    for _ in range(k): x = mul(x[0], x[1], a, b)
    return x
rows = defaultdict(list)
for i in range(len(P)):
    for k in range(i + 1, len(P)):
        P1, P2 = P[i], P[k]
        if P1.q == P2.q: continue
        ma, mb = mul(P1.a, P1.b, P2.a, P2.b); N = P1.N * P2.N
        # gamma_2(m) directly: residues x mod N represent O/m
        x = np.arange(N, dtype=np.int64)
        j1 = P1.symbol(x, np.zeros_like(x)); j2 = P2.symbol(x, np.zeros_like(x))
        ok = (j1 >= 0) & (j2 >= 0)
        ch2 = np.where(ok, Z6[(2 * (j1 + j2)) % 6], 0)
        g2m = np.sum(ch2 * np.exp(-2j * np.pi * ((mb * x) % N) / N)) / np.sqrt(N)
        al = cval(ma, mb); al = al / abs(al)
        W = np.conj(al) ** 2
        for Pr, other in ((P1, P2), (P2, P1)):
            # sigma_p = lambda^2 * other ; lambda^3 m/p = lambda^3 * other ; eps_p = -(lambda^5 other^2)^{-1}
            s_a, s_b = mul(*powel(*lam, 2), other.a, other.b)
            e_a, e_b = mul(*powel(*lam, 5), *mul(other.a, other.b, other.a, other.b))
            js = sym(Pr, s_a, s_b)
            # chi_p(eps)^{-3} = chi_p(-1)^{-3} chi_p(lambda^5 other^2)^{+3}
            je = (-3 * sym(Pr, -1, 0) + 3 * sym(Pr, e_a, e_b)) % 6
            W *= Z6[(-2 * js) % 6] * gam_p(Pr, 1, -1) * gam_p(Pr, 3, +1) * Z6[je]
        rows[(ma % 72, mb % 72)].append((W, W * g2m ** 4))
nc = bad = 0; spreadW = []
for key, L in rows.items():
    if len(L) < 2: continue
    nc += 1
    V = np.array([v for _, v in L]); Wv = np.array([w for w, _ in L])
    if np.max(np.abs(V - V[0])) > 1e-6: bad += 1
    spreadW.append(np.max(np.abs(Wv - Wv[0])))
print("composite m=p1p2: classes mod 72 with >=2 m: %d ; classes where W*gamma2(m)^4 not constant: %d ; median spread of W: %.3f"
      % (nc, bad, float(np.median(spreadW))))
print("total m:", sum(len(L) for L in rows.values()))
