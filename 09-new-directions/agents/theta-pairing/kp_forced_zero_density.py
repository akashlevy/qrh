"""Upper bound for the natural density of rows m (coprime to 6) at which the Kazhdan-Patterson
relations FORCE rho(m)=0:  v_p(m) = 5 mod 6  (tau(m p^5)=0), or v_p(m) = 2 mod 6 with G_3 != 1
(we count ALL such m, i.e. assume G_3 != 1 always -> upper bound).  Primes of Z[omega] prime to 6:
split p = 1 mod 3 -> two ideals of norm p; inert p = 2 mod 3, p != 2 -> one ideal of norm p^2."""
from sympy import primerange
def local_forced(N):
    # P(v = k) = (1-1/N) N^-k ; forced if k = 2 or 5 mod 6
    s = 0.0
    for k in range(2, 200):
        if k % 6 in (2, 5):
            s += (1 - 1 / N) * N ** (-k)
    return s
prod_free = 1.0
for p in primerange(5, 2_000_000):
    if p % 3 == 1:
        f = local_forced(p); prod_free *= (1 - f) ** 2
    elif p % 3 == 2:
        prod_free *= (1 - local_forced(p * p))
print("density of rows with a KP-forced zero  <=", 1 - prod_free)
print("CS gain factor from removing them      >=", prod_free ** 0.5)
