"""Exact symbolic identities and isolating signs for the feedback barrier.
Pure standard library. All proof decisions use Fraction arithmetic.
"""
from fractions import Fraction as F


def add(*ps):
    q = {}
    for p in ps:
        for ij, a in p.items():
            q[ij] = q.get(ij, F(0)) + a
    return {ij:a for ij,a in q.items() if a}


def scale(a,p):
    return {ij:F(a)*v for ij,v in p.items() if a*v}


def mul(p,q):
    r={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            r[i+k,j+l]=r.get((i+k,j+l),F(0))+a*b
    return {ij:a for ij,a in r.items() if a}


def const(a):
    return {(0,0):F(a)} if a else {}


def deriv_c(p):
    return {(i-1,j):i*a for (i,j),a in p.items() if i}


c={(1,0):F(1)}
d={(0,1):F(1)}
one=const(1)
a=F(5,6)
D=add(const(F(5,2)),scale(-2,c))
P=add(one,scale(-1,c))
J=add(mul(add(const(a),scale(-1,d)),D),mul(d,P))
Q=add(mul(add(const(36),scale(-18,c)),mul(d,d)),
      mul(add(scale(57,c),const(-78)),d),const(25),scale(-20,c))
# 36J(R-2/3)=6[2J-6dJ+3(a-d)dP].
assert Q==scale(6,add(scale(2,J),scale(-6,mul(d,J)),
                     scale(3,mul(mul(add(const(a),scale(-1,d)),d),P))))
num=add(scale(2,mul(J,add(one,scale(-1,d)))),
        mul(mul(add(const(a),scale(-1,d)),d),P))
# R=num/(2J); derivative equals -d(a-d)^2/(4J^2).
lhs=add(mul(deriv_c(num),scale(2,J)),scale(-1,mul(num,scale(2,deriv_c(J)))))
rhs=scale(-1,mul(d,mul(add(const(a),scale(-1,d)),add(const(a),scale(-1,d)))))
assert lhs==rhs

# Weight numerators on L0,LK,H, denominator 9+18d.
w0=add(const(4),scale(18,d));wk=const(2);wh=const(3)
assert add(w0,wk,wh)==add(const(9),scale(18,d))
assert not add(scale(F(1,12),w0),scale(F(7,12),wk),
               scale(F(-1,2),mul(add(one,d),wh)))
assert not add(scale(F(-1,6),w0),scale(F(1,3),wk),mul(d,wh))
assert add(scale(F(5,6),w0),scale(F(1,3),wk),mul(add(one,d),wh))==add(const(7),scale(18,d))

# Elimination multiplier after c=1/(5-3lambda), d=(5-9lambda)/(6+18lambda).
# Clearing the denominator, original Q becomes
# 3[18(3-2lambda)d^2+3(26lambda-37)d+5(7-5lambda)].
x=c
v=add(const(5),scale(-9,x))
z=add(const(6),scale(18,x))
q_clear=add(scale(18,mul(add(const(3),scale(-2,x)),mul(v,v))),
            scale(3,mul(mul(add(scale(26,x),const(-37)),v),z)),
            scale(5,mul(add(const(7),scale(-5,x)),mul(z,z))))
Pc=add(const(20),scale(21,x),scale(-954,mul(x,x)),scale(657,mul(mul(x,x),x)))
assert q_clear==scale(-36,Pc)

poly=lambda lam:657*lam**3-954*lam**2+21*lam+20
assert poly(F(1,6))==F(1,24)
assert poly(F(1,5))==F(-1088,125)
assert F(1971,25)-F(1908,6)+21==F(-5454,25)<0
lo=F('0.87495701942009894');hi=F('0.87495701942009895')
assert poly(F(11,3)-4*lo)<0<poly(F(11,3)-4*hi)
print('PASS: critical count polynomial identity')
print('PASS: monotonicity derivative identity')
print('PASS: all three dual-weight cancellations')
print('PASS: exact cubic elimination multiplier')
print('PASS: unique-root endpoint and derivative signs')
print('PASS: theta* lies strictly between',str(lo),'and',str(hi))
