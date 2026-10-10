# Referee check: the j = 0 Poisson node is a Ramanujan sum (-1 at unit H), i.e. it carries a Moebius sign
# with a q_p^{-1/2} deficit relative to the unimodular gamma_j(p) of j != 0; and its s-A cross term is
# conj(chi_r(p))^a (not cancellable unless a = 0).
import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
src = open("/Users/akashlevy/Documents/riemann/gram_split/vls/D_cross_term.py").read().split("P = [Q for Q in primes_upto(400)")[0]
exec(src)
P = [Q for Q in primes_upto(300) if Q.kind == "split"]
err = 0; n = 0
for a in [1, 2]:
    for R in P[:4]:
        for Pp in P[:6]:
            if R.N == Pp.N: continue
            for H in [(1, 0), (2, 1), (-1, 3)]:
                if sym1(R, *H) < 0 or sym1(Pp, *H) < 0: continue
                pred = R.N * Z6[sym1(R, *H)] ** a * np.conj(Z6[sym1(R, Pp.a, Pp.b)]) ** a * (-1.0)
                err = max(err, abs(Tbrute(R, Pp, a, 0, H) - pred) / abs(pred)); n += 1
print("j=0 node: T = q_r chi_r(H)^a conj(chi_r(p))^a * (-1)  [Ramanujan c_p(H) = -1]: %d cases, max rel err %.2e" % (n, err))
