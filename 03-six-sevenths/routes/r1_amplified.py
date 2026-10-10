"""
R1 variant with weights favouring sixth powers: prime amplifier inside the mean square.
  S_amp = sum_{0<N u<=Z^m} |M_u(Z^r) Q_u|^2,  M_u = Z^{-r/2} A_u(Z^r),  Q_u = Z^{-z/2} sum_{p ~ Z^z} c_p chi_p(u)
(paper78 Lemma lem:marked with one slot; c_p = 1).  Normalise m = 1.
Rows u = q^6 (q prime, N q ~ Y = Z^{1/6}) have |Q_u|^2 ~ Z^{z}, so if S_amp << Z^{1+eps}:
  |A_1(D)|^2 << D Z / (Y P) + D^2 / Y^2,  D = Z^r, P = Z^z
  => theta(r,z) = max( (r + 5/6 - z)/r , (2r - 1/3)/r ) / 2.
Constraint sets (all with r > 0, z >= 0):
  (L)   Lemma marked as proved:       r + 2z <= 1,  2r + 8z <= 3
  (W)   prior-analysis capacity under W(w): z <= min((1 + w/2 - r)/2, 1 - r, (3 - 2r)/8)
  (D8)  case-B wall removed:           r + z <= 1,   2r + 8z <= 3
  (D)   canonical diagonal only:       r + z <= 1     (true under GRH; 'lossless' marked moment)
  all sets additionally impose z <= r (P <= D).
"""
from fractions import Fraction as F
import itertools

def theta(r, z):
    return max((r + F(5, 6) - z) / r, (2 * r - F(1, 3)) / r) / 2

def feasible(r, z, kind, w=F(0)):
    if r <= 0 or z < 0 or z > r: return False      # slot no longer than the inverse sum (P <= D)
    if kind == "L":  return r + 2 * z <= 1 and 2 * r + 8 * z <= 3
    if kind == "W":  return z <= min((1 + w / 2 - r) / 2, 1 - r, (3 - 2 * r) / 8)
    if kind == "D8": return r + z <= 1 and 2 * r + 8 * z <= 3
    if kind == "D":  return r + z <= 1
    raise ValueError

def best(kind, w=F(0), N=240):
    bestv, arg = None, None
    for i in range(1, 2 * N + 1):
        r = F(i, N)
        for j in range(0, N + 1):
            z = F(j, N)
            if feasible(r, z, kind, w):
                t = theta(r, z)
                if bestv is None or t < bestv: bestv, arg = t, (r, z)
    return bestv, arg

for kind in ("L", "D8", "D"):
    v, (r, z) = best(kind)
    print(f"({kind:2s}) min theta on grid = {v} = {float(v):.6f} at r={r}, z={z}  (P = D^{float(z/r) if r else 0:.4f})")
for w in (F(1, 10), F(1, 5), F(1, 2)):
    v, (r, z) = best("W", w)
    print(f"(W) w={w}: min theta = {float(v):.6f} at r={r}, z={z}")
# exact checks
print("exact: (L) at r=1,z=0:", theta(F(1), F(0)), "; (D8) at r=5/6,z=1/6:", theta(F(5,6), F(1,6)),
      "; (D) at r=7/12,z=5/12:", theta(F(7,12), F(5,12)), "; (D) at r=1/2,z=1/2:", theta(F(1,2), F(1,2)))
# does the 6/7 point satisfy Lemma marked?  r + 2z and 2r + 8z at (7/12, 5/12):
r, z = F(7, 12), F(5, 12)
print(f"6/7 point (r,z)=(7/12,5/12): r+z={r+z}, r+2z={r+2*z} (needs <=1), 2r+8z={2*r+8*z} (needs <=3)")
# along the (L) boundary r + 2z = 1 the extraction exponent is 3/4 + 1/(6r) >= 11/12 (decreasing in r, r <= 1)
for r in (F(1, 2), F(2, 3), F(5, 6), F(1)):
    z = (1 - r) / 2
    assert theta(r, z) == F(3, 4) + 1 / (6 * r) or theta(r, z) == 1 - 1 / (6 * r)
    print(f"  (L) boundary r={r}, z={(1-r)/2}: theta={theta(r, z)}")
