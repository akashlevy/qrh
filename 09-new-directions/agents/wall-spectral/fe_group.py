"""Check the functional-equation group of the energy MDS (NOTES §16): finite or infinite?"""
import numpy as np, itertools
s1 = np.array([[-1,0,0],[0,1,0],[2,0,1]])
s2 = np.array([[1,0,0],[0,-1,0],[0,2,1]])
s3 = np.array([[1,0,1],[0,1,1],[0,0,-1]])
gens = [s1,s2,s3]
for g in gens: assert (g@g == np.eye(3)).all()
seen = {tuple(np.eye(3,dtype=int).ravel())}
frontier = [np.eye(3,dtype=int)]
for depth in range(1,13):
    new = []
    for M in frontier:
        for g in gens:
            P = g@M
            t = tuple(P.ravel())
            if t not in seen:
                seen.add(t); new.append(P)
    frontier = new
    print(depth, len(seen), max(abs(np.array(t)).max() for t in seen))
# invariant quadratic form?
import sympy as sp
a,b,c,d,e,f = sp.symbols('a b c d e f')
Q = sp.Matrix([[a,d,e],[d,b,f],[e,f,c]])
eqs = []
for g in gens:
    G = sp.Matrix(g)
    eqs += list(G.T*Q*G - Q)
print(sp.solve(eqs,[a,b,c,d,e,f],dict=True))
# fixed vectors / invariant functionals
for name,g in zip("123",gens):
    print(name, np.linalg.eigvals(g))
w = s1@s2@s3
print("Coxeter element eigenvalues", np.linalg.eigvals(w))
