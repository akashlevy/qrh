"""Cubic architecture in an enlarged geometry box (Q >> Z^2 allowed), conservative inputs:
hmin = 3/2 (cubic P2 mean square from the cubic large sieve alone, no Rankin-Selberg/Kubota FE),
cB = 4 (sextic Kubota index h^4, rigorous upper bound) -- the winning node is the trivial Poisson node."""
import sys, json, numpy as np
import model_cubic as MC
def solve_box(cfg_over, lx_max, ly_max, Smax, nS=20):
    orig = MC.linprog
    def lp2(*a, **kw):
        kw["bounds"] = [(0, lx_max), (0, ly_max), (None, None)]; return orig(*a, **kw)
    MC.linprog = lp2
    try:
        return MC.theta(cfg_over, Sgrid=np.linspace(0.02, Smax, nS))
    finally:
        MC.linprog = orig
if __name__ == "__main__":
    over = json.loads(sys.argv[1]); lxm, lym, Smax = map(float, sys.argv[2:5])
    sg, g = solve_box(over, lxm, lym, Smax)
    print(json.dumps(dict(cfg=over, box=[lxm, lym, Smax], theta=sg,
                          geom={k: (float(v) if not isinstance(v, tuple) else v) for k, v in g.items()})))
