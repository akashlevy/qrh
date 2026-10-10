"""Exact identities for explicitly delimited QRH exponent-template barriers.

No claim is made about possible stronger estimates or actual zero locations.
"""
from fractions import Fraction as F
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "06-kappa-bootstrap", "certificates"))  # repo layout
from bootstrap_certificate import add, mul, scale, constant


def scaled_sum(*terms):
    return [sum(scale*v[i] for scale, v in terms) for i in range(len(terms[0][1]))]


# Variables M, ell, e, constant. With the unchanged Gram majorant,
# theta >= 5/6 + M/12 - ell/6 + e/2. With the two retained energy losses,
# e >= 0, e >= (1+3ell-2M)/4, e >= M+ell-1.
A = [F(1, 2), -F(3, 4), F(1), -F(1, 4)]
B = [-F(1), -F(1), F(1), F(1)]
E = [F(0), F(0), F(1), F(0)]
assert scaled_sum((F(1, 5), A), (F(1, 60), B), (F(17, 60), E)) == [
    F(1, 12), -F(1, 6), F(1, 2), -F(1, 30)]
assert F(5, 6)+F(1, 30) == F(13, 15)
print('PASS: unchanged low-majorant template requires theta >= 13/15.')

# Prime-sixth-power extraction: theta=max(1/2+5k/12,1-k/6).
# A convex combination cancels k and equals 6/7.
assert F(2, 7)*F(5, 12)+F(5, 7)*(-F(1, 6)) == 0
assert F(2, 7)*F(1, 2)+F(5, 7) == F(6, 7)
assert F(1, 2)+F(5, 12)*F(6, 7) == 1-F(6, 7)/6 == F(6, 7)
print('PASS: prime-sixth-power extraction has exact minimum theta=6/7.')


def feedback_poly(ell):
    return 657*ell**3-954*ell**2+21*ell+20


# Independently verify the cleared substitution in the critical-class bound.
# c=1/H, delta=N/D, ell is the polynomial variable.
l = {(1, 0): F(1)}
H = add(constant(5), scale(-3, l))
N = add(constant(5), scale(-9, l))
D = add(constant(6), scale(18, l))
Q_cleared = add(mul(add(scale(36, H), constant(-18)), mul(N, N)),
                mul(add(constant(57), scale(-78, H)), mul(N, D)),
                mul(add(scale(25, H), constant(-20)), mul(D, D)))
P = {(3, 0): F(657), (2, 0): F(-954), (1, 0): F(21), (0, 0): F(20)}
assert Q_cleared == scale(-108, P)
print('PASS: exact critical-class substitution Q=-3P/[H(1+3ell)^2].')

# Three positive weights before normalization by 9+18delta.
# Coefficients in the basis 1,delta; combinations eliminate M and ell.
def affine_sum(*terms):
    return tuple(sum(coef*v[i] for coef, v in terms) for i in (0, 1))
w0, wK, wH = (F(4), F(18)), (F(2), F(0)), (F(3), F(0))
assert affine_sum((1, w0), (1, wK), (1, wH)) == (9, 18)
assert affine_sum((F(1,12), w0), (F(7,12), wK), (-F(3,2), (1,1))) == (0,0)
assert affine_sum((-F(1,6), w0), (F(1,3), wK), (F(3), (0,1))) == (0,0)
assert affine_sum((F(5,6), w0), (F(1,3), wK), (F(3), (1,1))) == (7,18)
print('PASS: positive dual weights yield theta >= (7+18delta)/(9+18delta).')

# A critical-class necessary condition for ideal-capacity feedback.
# The dual argument itself allows unbalanced geometry. Isolation alone is
# not a proof that its boundary is attained by an analytic argument.
left, right = F(16683858898627, 10**14), F(16683858898628, 10**14)
assert feedback_poly(left) > 0 > feedback_poly(right)
# Derivative negative on [1/6,1/5] using an elementary uniform bound.
assert 1971*F(1, 5)**2-1908*F(1, 6)+21 < 0
print('PASS: exactly one feedback root in [1/6,1/5], isolated at', left, right)
print('Corresponding theta lies in', F(11, 12)-right/4, F(11, 12)-left/4)
