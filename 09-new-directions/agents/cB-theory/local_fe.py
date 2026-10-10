"""
Local functional equation of the COPRIME parts of the n-th order Kubota series at a prime p | index.
(cB-theory, 10 Oct 2026.)  Exact symbolic algebra (sympy); the only analytic input is the rank-one
Kubota FE of the FULL series  psi(r,u) = N(r)^{1/2-u} * G(u) * [Z(1-u)/Z(u)] * psi(r,1-u)
(vector-valued over S-classes, S-matrix independent of r except through r mod F_S^{x n}),
with Z(u) = zeta_F(n u - n/2 + 1) and u the normalised variable (g~ = g/N^{1/2}).

Notation at one prime p (q = N(p), X = q^{-u}, Y = q^{-(1-u)} = 1/(q X)):
  C_j(u) = sum_{(c,p)=1} (p/c)^{-j} g~(m',c) N(c)^{-u}   (j mod n)  = coprime part at index p^j m'.
  psi_l  = full series at index p^l m'.
p-parts (g~(p^l,p^k) for k<=l, n|k : q^{k/2}(1-1/q);  k = l+1 : q^{l/2} G_{l+1}, or -q^{(l-1)/2} if n|l+1;
all other k : 0) and twisted multiplicativity give, for 0 <= l <= n-1,
  psi_l = C_l + a_l X^{l+1} C_{l'},   l' = -l-2 mod n,   |a_l| = q^{l/2}   (a_{n-1} = -q^{(n-2)/2}).
We solve for the FE  T(u) = G(u) M(u) T(1-u)  of T_j = Z^{(p)} C_j  and print every monomial of M
as  weight * q^{e(1/2-u)}  (weight = size on Re u = 1/2, e = local conductor exponent).
"""
import sympy as sp

q, X = sp.symbols('q X', positive=True)
Y = 1 / (q * X)


def monomials(expr):
    """expr is a Laurent polynomial in X with coefficients in q^(rational) * (unit symbols)."""
    expr = sp.expand(expr)
    out = []
    for term in sp.Add.make_args(expr):
        c, rest = term.as_independent(X)
        b = sp.simplify(X * sp.diff(rest, X) / rest) if rest.has(X) else 0
        c = sp.simplify(term / X ** b)
        # weight on Re u = 1/2:  c * q^{-b/2};  conductor exponent b
        w = sp.simplify(c * q ** sp.Rational(-b, 2))
        out.append((b, w))
    return sorted(out, key=lambda t: t[0])


def block(n, l, verbose=True):
    lp = (-l - 2) % n
    D = lambda Z: 1 - q ** sp.Rational(n - 2, 2) * Z ** n
    if lp == l:   # singleton
        if l == n - 1:
            b = -q ** sp.Rational(n - 2, 2)
        else:
            e = sp.Symbol(f'e{l}')            # unit phase (here a ray/quadratic Gauss-sum value)
            b = q ** sp.Rational(l, 2) * e
        A = lambda Z: 1 + b * Z ** (l + 1)
        # T_l(u) = q^{l/2} X^l * A(Y)/A(X) * D(X)/D(Y) * T_l(1-u)   [Z^{(p)} normalisation]
        num = q ** sp.Rational(l, 2) * X ** l * A(Y)
        if verbose:
            print(f"n={n} l={l}: SINGLETON.  psi_l = (1 + b X^{l+1}) C_l,  |b| = q^{sp.Rational(l,2) if l<n-1 else sp.Rational(n-2,2)}")
            print(f"   M = q^(l(1/2-u)) * [A(Y)/A(X)] * [D(X)/D(Y)]  with A,D Euler-type polynomials")
            print(f"   leading monomials of q^(l/2) X^l A(Y): {monomials(num)}")
        return {('d', l): monomials(num)}
    # pair
    el, elp = sp.symbols(f'e{l} e{lp}')
    al = q ** sp.Rational(l, 2) * el
    alp = q ** sp.Rational(lp, 2) * elp
    A = lambda Z: sp.Matrix([[1, al * Z ** (l + 1)], [alp * Z ** (lp + 1), 1]])
    F = sp.diag(q ** sp.Rational(l, 2) * X ** l, q ** sp.Rational(lp, 2) * X ** lp)
    # T(u) = adj(A(X)) F A(Y) / D(Y) * T(1-u)   (det A(Z) = D(Z) once e_l e_l' = 1)
    M = (A(X).adjugate() * F * A(Y)).subs(elp, 1 / el)
    M = M.applyfunc(sp.expand)
    detchk = sp.simplify(A(X).det().subs(elp, 1 / el) - D(X))
    res = {}
    if verbose:
        print(f"n={n} pair (l,l')=({l},{lp}):  det A - D = {detchk}")
    for i, a in enumerate((l, lp)):
        for j, b in enumerate((l, lp)):
            mons = monomials(M[i, j])
            res[(a, b)] = mons
            if verbose:
                print(f"   T_{a} <- T_{b}(1-u):  " + ",  ".join(f"[e={e}, weight {w}]" for e, w in mons))
    return res


if __name__ == "__main__":
    for n in (2, 3, 6):
        print("=" * 100)
        for l in range(n):
            block(n, l)
    print("=" * 100)
    print("Reading: the coefficient of T_l(1-u) in the FE of T_l is q^{l(1/2-u)}(1 - q^-1) (weight 1,")
    print("exponent l); the only other monomials are off-diagonal, of weight q^{-1/2}, exponents -1 and n-1,")
    print("i.e. the shifted series a_l X^{l+1} psi_{l'} of conductor q^{l'}.  1/D(Y) only adds exponents -n k")
    print("with weights q^-k.  Hence: effective conductor of the coprime part at index p^l is q^{(l mod n)}.")
