"""theta with mark amplitude x capped at xmax (diag copy of coord/xcap.py, with caseA/delta0 overrides)."""
import sys, json, numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/diag")
import model as M
_orig = np.linspace
def run(xmax, over):
    def lin(a, b, n, *args, **kw):
        if a == 0 and b == 0.5 and n == M.DEFAULT["xgrid"]:
            return _orig(0, xmax, n)
        return _orig(a, b, n, *args, **kw)
    M.np.linspace = lin
    try:
        return M.theta(over)
    finally:
        M.np.linspace = _orig
over = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
sg, g = run(float(sys.argv[1]), over)
print(json.dumps(dict(xmax=float(sys.argv[1]), cfg=over, theta=sg, geom=g)))
