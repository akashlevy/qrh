import numpy as np
exec(open("k0.py").read().split("# ---------- C sum")[0])
for (mu, qd, dch) in dlist:
    e = Etilde(0 * 4 + 1 if False else (1 % 4) * 4 + 0, dch)
    print("d with q_d=%2d  mu=%+d  Etilde(C=1,d) = %.6f%+.6fi   weighted mu q_d^{-5/6} Etilde = %.6f" % (qd, mu, e.real, e.imag, (mu * qd ** (-5 / 6) * e).real))
tot = sum(mu * qd ** (-5 / 6) * Etilde(4, dch) for (mu, qd, dch) in dlist)
print("sum_d (C=1) =", tot)
print("IOm = int Omega(u) u^{-1/6} du =", IOm, "  IH =", IH, "  ghat0 =", ghat0)
# E without the d-twist and without units collapse, for orientation: is the C=1 ray average positive?
