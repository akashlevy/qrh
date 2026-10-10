"""gamma_1(s) from probe.svec (mu G conj beta) vs brute-force sextic Gauss sums over Z[w]/s."""
import math, numpy as np
from nlib import Arith, W
from probe import svec
ar = Arith(410000)
ar.k4 = np.array([int(p.sextic_exp(4, 0)[0]) for p in ar.P])
def brute(fac, A, B):
    Nn = A * A - A * B + B * B
    a = np.arange(Nn); X, Y = np.meshgrid(a, a, indexing="ij"); X = X.ravel(); Y = Y.ravel()
    e = np.zeros(len(X), np.int64); z = np.zeros(len(X), bool)
    for i in fac:
        k = ar.P[i].sextic_exp(X, Y).astype(np.int64); z |= k < 0; e += np.where(k < 0, 0, k)
    zz = (X + Y * W) / complex(A + B * W)
    v = np.where(z, 0, np.exp(2j * np.pi * (e % 6) / 6)) * np.exp(4j * np.pi * zz.imag / math.sqrt(3))
    return v.sum() / Nn / math.sqrt(Nn)
worst = 0; n = 0
for Y in (128, 256, 512):
    A, B, Nn, fac, coef, cl = svec(ar, Y)
    w = np.sin(2 * np.pi * (Nn / Y - 0.5)) ** 2
    for j in range(0, len(A), max(1, len(A) // 25)):
        if w[j] < 1e-3: continue
        g = brute(fac[j], int(A[j]), int(B[j]))
        worst = max(worst, abs(coef[j] / w[j] - g)); n += 1
print(f"gamma_1 identity vs brute force on {n} squarefree s (N <= 512, incl. composites/inert): max |diff| = {worst:.2e}")
