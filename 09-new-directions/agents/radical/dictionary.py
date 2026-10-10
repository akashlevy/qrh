"""radical: exact dictionary between a zero-free half-plane Re s > theta and the
equivalent (or one-way) statements of the 'radical' reformulations, evaluated at
theta = 7/8, the pending record 0.874887, 27/31 and 6/7.  Also the bootstrap
(H(beta*)) Halasz-Montgomery threshold.  Pure arithmetic; see RADICAL.md for proofs/refs.
"""
from fractions import Fraction as Fr

REC = 683505193 / 781250000          # tracker PR #9 (pending)
cases = [("7/8", 7 / 8), ("record 0.874887", REC), ("27/31", 27 / 31), ("6/7", 6 / 7), ("19/22", 19 / 22)]

print(f"{'theta':>18} | {'horocycle a':>11} | {'Beurling p':>10} | {'res. gap g':>10} | {'LP decay':>8} | {'dBN Lambda<=':>12} | {'HM boot thr':>11}")
for name, th in cases:
    a = 1 - th / 2                   # Zagier: a_0(y) - (3/pi) int f = O(y^{a-eps}) for all f  <=>  no zeros in Re s > 2(1-a)
    p = 1 / th                       # Beurling 1955: B_p = L^p(0,1)  <=>  zeta != 0 on Re s > 1/p
    g = (1 - th) / 2                 # SL2(Z) resonances at rho/2: strip-free width 1/2 - beta/2
    lp = g                           # Lax-Phillips semigroup decay rate (spectral bound) = g
    lam = 2 * (th - 0.5) ** 2        # de Bruijn: |Im z| <= 2 theta - 1  =>  Lambda <= (2theta-1)^2/2  (ONE-WAY)
    hm = (1 + th) / 2                # Halasz-Montgomery under H(beta*=theta): effective only for sigma > (1+beta*)/2
    print(f"{name:>18} | {a:11.5f} | {p:10.5f} | {g:10.5f} | {lp:8.5f} | {lam:12.5f} | {hm:11.5f}")

print()
print("Unconditional endpoints before P78 (all = the trivial theta = 1):")
print("  horocycle a = 1/2 (Sarnak 1981, continuous spectrum on Re s = 1/2); Beurling p = 1 (Bercovici-Foias);")
print("  resonance gap 0; dBN Lambda <= 1/2 (de Bruijn), < 1/2 (Ki-Kim-Lee), 0.2 (Platt-Trudgian, numerical).")
print("To beat Platt-Trudgian Lambda <= 0.2 through the one-way conversion one needs theta <= %.5f" % (0.5 + (0.1) ** 0.5))
print("To beat Polymath15 Lambda <= 0.22 one needs theta <= %.5f" % (0.5 + (0.11) ** 0.5))

# Halasz-Montgomery large-values count under H(beta*): R V^2 <= N (N + R C N^{beta*+eps} Q^eps), V = N^sigma
# -> R << N^{2-2 sigma} as soon as N^{1+beta*-2 sigma} < 1/2, i.e. sigma > (1+beta*)/2.
# Unconditional off-diagonal (smoothed sum via convexity on Re s = 1/2): N^{1/2} Q^{1/4}:
# -> needs N^{3/2-2 sigma} Q^{1/4} < 1, impossible for sigma < 3/4.
print()
for bs in (0.875, REC, 27 / 31):
    print(f"H(beta*={bs:.6f}): HM count effective only for sigma > {(1 + bs) / 2:.5f};"
          f" P78's binding rows sit at a = (1+delta)/2 in [0.63, 0.70]  -> no effect")
