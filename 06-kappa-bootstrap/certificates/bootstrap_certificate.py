"""Exact arithmetic for a proposed kappa-bootstrap of the QRH argument.

This proves the stated polynomial and scale inequalities, NOT the analytic
moment theorem. The accompanying proof must justify its extended kappa range.
No floating-point number is used in a proof decision. Python standard library.
Run with --negative-control to require rejection of an overstrong target.
"""
from fractions import Fraction as F
from math import comb
from itertools import product
import argparse
import json


def add(*polys):
    out = {}
    for p in polys:
        for ij, a in p.items():
            out[ij] = out.get(ij, F(0)) + a
    return {ij: a for ij, a in out.items() if a}


def scale(a, p):
    return {ij: a * b for ij, b in p.items() if a * b}


def mul(p, q):
    out = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            ij = i + k, j + l
            out[ij] = out.get(ij, F(0)) + a * b
    return {ij: a for ij, a in out.items() if a}


def constant(a):
    return {(0, 0): F(a)}


def evaluate(p, d, x):
    return sum(a * d**i * x**j for (i, j), a in p.items())


def bernstein_bounds(p, box):
    """Convex hull bounds on a rectangle, after exact affine substitution."""
    d0, d1, x0, x1 = box
    n = max(i for i, j in p)
    m = max(j for i, j in p)
    power = {}
    for (i, j), a in p.items():
        for k in range(i + 1):
            for l in range(j + 1):
                value = (a * comb(i, k) * d0**(i-k) * (d1-d0)**k
                         * comb(j, l) * x0**(j-l) * (x1-x0)**l)
                power[k, l] = power.get((k, l), F(0)) + value
    values = [sum(power.get((k, l), F(0))
                  * F(comb(i, k), comb(n, k)) * F(comb(j, l), comb(m, l))
                  for k in range(i+1) for l in range(j+1))
              for i, j in product(range(n+1), range(m+1))]
    return min(values), max(values)


def certify(p, box, depth=0, leaves=None):
    if leaves is None:
        leaves = []
    lower, upper = bernstein_bounds(p, box)
    if lower >= 0:
        leaves.append((box, lower))
        return leaves
    assert depth < 46, ('Undecided positivity', box, lower, upper)
    assert upper >= 0, ('Negative rectangle', box, upper)
    d0, d1, x0, x1 = box
    # First reject actual negative points, avoiding huge subdivision on failure.
    for d, x in product((d0, (d0+d1)/2, d1), (x0, (x0+x1)/2, x1)):
        assert evaluate(p, d, x) >= 0, ('Negative value', d, x)
    if d1-d0 >= x1-x0:
        mid = (d0+d1)/2
        halves = [(d0, mid, x0, x1), (mid, d1, x0, x1)]
    else:
        mid = (x0+x1)/2
        halves = [(d0, d1, x0, mid), (d0, d1, mid, x1)]
    for half in halves:
        certify(p, half, depth+1, leaves)
    return leaves


def build_certificate(theta, kappa, b, margin):
    ell = 4 * (F(11, 12) - theta)
    M = 1 - ell
    lx, ly = (M-b)/2, (M+b)/2
    h = 1-lx+ell
    c = 1/(6*kappa)
    alpha, z0, d0 = F(5, 6), F(17, 50), F(1, 50)
    delta, x = {(1, 0): F(1)}, {(0, 1): F(1)}
    Dx = add(constant(3), scale(-(1+4*c), x))
    Px = mul(add(constant(1), scale(-1, x)), add(constant(2), scale(-4*c, x)))
    J = add(mul(add(constant(alpha), scale(-1, delta)), Dx), mul(delta, Px))
    # E=(1+d)(1-ly)/2-ell/2+x*d*ell+h*(R+d/2-1/6)-theta.
    # R=1-d+(alpha-d)*d*P/(2J). N=J*(-E-margin).
    bracket = add(constant(theta-(1-ly)/2+ell/2-h*F(5, 6)-margin),
                  scale(-(1-ly)/2+h/2, delta), scale(-ell, mul(delta, x)))
    N = add(mul(J, bracket),
            scale(-h/2, mul(mul(add(constant(alpha), scale(-1, delta)), delta), Px)))
    box = (F(0), alpha, F(0), F(1, 2))
    checks = {}

    def require(name, value, strict=True):
        assert value > 0 if strict else value >= 0, (name, value)
        checks[name] = str(value)

    require('kappa >= 13/18', kappa-F(13, 18), False)
    require('kappa <= 3/4', F(3, 4)-kappa, False)
    require('positive b', b)
    require('positive ell', ell)
    require('X rescaled length', lx-ell)
    require('Y rescaled length', ly-ell)
    require('Q rescaled length', M-2*ell)
    require('Gram domination', ly-ell-F(11, 6)*b)
    require('unrescaled energy A', 2-10*ell)
    # f(d)=-d+max(0,5*ell-1+d)/8 is convex, so endpoints suffice.
    for d in (F(0), ell):
        require('rescaled energy tuple d='+str(d),
                d-max(F(0), 5*ell-1+d)/8, False)
    require('prime supply', ell/h-F(7, 37))
    require('row extension supply', 5*ell-h)
    require('small rows', -(h*(z0-F(1, 6))-ly/2+F(63, 5000)))
    require('w residue margin', ly/20)
    require('z residue margin', h/600)
    # At sigma >= theta, w>=19/20,z>=33/200.
    z, w = F(33, 200), F(19, 20)
    for i, exponent in enumerate((-theta, -6*z, 4-5*theta-6*z, 1-w-6*z)):
        require('local error '+str(i), -theta-exponent, False)
    require('Euler defect A', -1-(4-6*theta-6*z))
    require('Euler defect B', -1-(1-theta-w-6*z))
    require('outer Euler region', theta+F(1, 2)-1)

    def E(d, q, R, row):
        a = (1+d)/2
        return (a-theta+h*(z0-F(1, 6))-a*ly-(1-a)*ell
                -(d/2-q)*ell+row*(R+d/2-z0))

    require('floor bin', -E(d0, d0/2, F(1), h))
    for d in (d0, alpha):
        require('intermediate rows '+str(d),
                -E(d, d/2, F(76, 75)-2*d/3, F(1, 2)))
    require('intermediate positive slope', F(101, 150)-alpha/6)
    require('J positive on full rectangle', bernstein_bounds(J, box)[0])
    # Independent direct checks of the polynomial transcription, not proof sampling.
    for d, xv in product((F(0), F(1, 3), F(2, 5), alpha), (F(0), F(1, 4), F(1, 2))):
        j = (alpha-d)*(3-(1+4*c)*xv)+d*(2-4*c*xv)*(1-xv)
        r = 1-d+(alpha-d)*d*(2-4*c*xv)*(1-xv)/(2*j)
        assert evaluate(N, d, xv) == j*(-E(d, xv*d, r, h)-margin)
    leaves = certify(N, box)
    area = sum((b1-b0)*(x1-x0) for (b0, b1, x0, x1), _ in leaves)
    assert area == alpha/2
    for subbox, lower in leaves:
        assert bernstein_bounds(N, subbox)[0] == lower >= 0
    return dict(theta=str(theta), decimal_theta=float(theta), kappa=str(kappa),
                c=str(c), b=str(b), ell=str(ell), lx=str(lx), ly=str(ly), h=str(h),
                endpoint_margin=str(margin), rectangle=[str(v) for v in box],
                checks=checks, bernstein_leaf_count=len(leaves),
                leaves=[dict(box=[str(v) for v in subbox], lower=str(lo))
                        for subbox, lo in leaves])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    parser.add_argument('--stage', choices=('1', '2', 'all'), default='all')
    parser.add_argument('--negative-control', action='store_true')
    args = parser.parse_args()
    # Stage 1 needs only the rational bound 34999/40000, not Liu's algebraic one.
    theta = F(17499141, 20000000)  # 0.87495705
    kappa = F(14999, 20000)       # 2*(34999/40000)-1
    if args.negative_control:
        # Below the critical class at this same kappa, hence must fail.
        try:
            build_certificate(F(8749569, 10000000), kappa, F(6170279, 50000000), F(1, 1000000000))
        except AssertionError as e:
            print('PASS: overstrong target rejected:', e)
        else:
            raise AssertionError('Negative control unexpectedly accepted')
    else:
        stages = [
            (theta, kappa, F(6170279, 50000000), F(1, 100000000)),
            (F(43747851, 50000000), F(7499141, 10000000),
             F(24681137, 200000000), F(1, 1000000000)),
        ]
        prior_bound = F(34999, 40000)
        results = []
        for number, parameters in enumerate(stages, 1):
            target, k, imbalance, reserve = parameters
            assert k == 2*prior_bound-1 and target < prior_bound
            if args.stage in ('all', str(number)):
                entry = build_certificate(*parameters)
                entry['stage'] = number
                entry['prior_bound_required'] = str(prior_bound)
                results.append(entry)
            prior_bound = target
        result = results if args.stage == 'all' else results[0]
        if args.output:
            with open(args.output, 'w') as out:
                json.dump(result, out, indent=2)
        print(json.dumps([{k: v for k, v in entry.items() if k not in ('leaves', 'checks')}
                          for entry in results], indent=2))
