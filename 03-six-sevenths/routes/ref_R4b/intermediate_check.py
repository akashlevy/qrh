# Referee: the paper's own intermediate-row bound (d_min <= d <= 1/2, t=1, no selected primes,
# R = 76/75 - 2delta/3, q = delta/2), eq. adaptive-intermediate, evaluated at general geometry.
from fractions import Fraction as F
def E(lx, ly, ell, s0, d, delta):
    h = 1 - lx + ell; a = (1 + delta) / 2; q = delta / 2
    R = F(76, 75) - F(2, 3) * delta
    return a * (1 - ly) - s0 - ell / 2 + q * ell + h * (F(17, 50) - F(1, 6)) + d * (R + delta / 2 - F(17, 50))
geoms = {"W(2/7)+refined": (F(11489, 28328), F(14840, 27149), F(32673, 97961)),
         "vertex A (4/7,4/7,3/7)": (F(4, 7), F(4, 7), F(3, 7)),
         "alt W(4/7) (22/49,34/49,3/7)": (F(22, 49), F(34, 49), F(3, 7)),
         "paper (17/48,23/48,1/6) at 7/8": (F(17, 48), F(23, 48), F(1, 6))}
for name, (lx, ly, ell) in geoms.items():
    s0 = F(7, 8) if "paper" in name else F(6, 7)
    vals = [E(lx, ly, ell, s0, F(1, 2), dl) for dl in (F(1, 50), F(3, 4), F(5, 6))]
    print(name, "E(1/2) at delta=1/50,3/4,5/6:", [float(v) for v in vals])
