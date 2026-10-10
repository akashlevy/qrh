"""N1.  Dunn-Radziwill Lemma 5.2 twist transform for the cubic theta, with psi = chi_r^j (r a split prime).
Paper conventions (zomega): chi_r(x) = (x/r)_6, cubic symbol (x/r)_3 = chi_r(x)^2, e(x/r) as in gamma_prime.
DR's Kloosterman transform (their (5.22)-(5.23)), written with the paper's additive character:
    K_j(u) = sum_{d mod r, (d,r)=1} psihat(d) (a/r)_3 e(u a / r),   a = (lambda^8 d)^{-1} mod r,
    psihat(d) = sum_x psi(x) e(x d / r).
Claim [R, algebra]:  K_j(u) = chi_r(lambda^8)^j * tau_j * tau_{j+2} * chi_r(u)^{-(j+2)}   (j, j+2 not = 0 mod 6),
    with tau_k = sum_x chi_r(x)^k e(x/r);  for j = 4 the a-sum is a Ramanujan sum; for j = 0, psihat(d) = -1.
So the reflected row picks up gamma_j gamma_{j+2} and the dual column factor chi^{-(j+2)}:
exactly the paper's omega_{p,j} = tau^-_{p,j} tau^+_{p,j+2} (..) and B_p(x) = chi_p(x)^{-j-2} (eq. reflection-row-phase)."""
import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *

def run(P):
    p = P.q
    x = np.arange(1, p, dtype=np.int64)
    sx = P.symbol(x, np.zeros_like(x))                      # exponent of chi(x), x = integer residue
    e = lambda t: np.exp(-2j * np.pi * ((P.b * (t % p)) % p) / p)   # e(t/r) for integer t (paper convention)
    chi = lambda k: Z6[(k * sx) % 6]                         # chi(x)^k on x = 1..p-1
    tau = {k: np.sum(chi(k) * e(x)) for k in range(6)}
    # lambda = 1 + 2 omega  ->  integer residue 1 + 2 r mod p
    lam = (1 + 2 * P.r) % p
    lam8 = pow(lam, 8, p)
    s_lam8 = int(P.symbol(np.array([lam8]), np.array([0]))[0])
    worst = 0.0
    for j in range(6):
        psihat = np.array([np.sum(chi(j) * e(x * d)) for d in x])        # d = 1..p-1
        a = np.array([pow(int(lam8 * d % p), -1, p) for d in x])
        sa = P.symbol(a, np.zeros_like(a))
        for u in [uu for uu in (1, 2, 5, 7, 11, p - 3) if uu % p]:
            K = np.sum(psihat * Z6[(2 * sa) % 6] * e(u * a))
            su = int(P.symbol(np.array([u % p]), np.array([0]))[0])
            if (j + 2) % 6 != 0:
                pred = Z6[(j * s_lam8) % 6] * (tau[j] if j % 6 else -1) * tau[(j + 2) % 6] * Z6[(-(j + 2) * su) % 6]
            else:   # j = 4: Ramanujan sum c_r(u) = -1 (u coprime)
                pred = Z6[(j * s_lam8) % 6] * tau[j] * (-1)
            worst = max(worst, abs(K - pred) / np.sqrt(p) ** 2)
    return worst

if __name__ == "__main__":
    Ps = [P for P in primes_upto(400) if P.kind == "split"][:12]
    print("max relative error of K_j(u) vs chi(lam^8)^j tau_j tau_{j+2} chi(u)^{-(j+2)} over j=0..5, 12 split primes:",
          max(run(P) for P in Ps))
    # the row factor gamma_j gamma_{j+2} for each j: which are ray (|.|=1 and constant on classes) ?
    print("j : row Gauss factor gamma_j*gamma_{j+2} ; dual column factor chi^{-(j+2)}")
    for j in range(6):
        print(j, ": gamma_%d gamma_%d ; chi^%d" % (j, (j + 2) % 6, (-(j + 2)) % 6))
