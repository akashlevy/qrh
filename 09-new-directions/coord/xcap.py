import sys, json, numpy as np
sys.path.insert(0, "/home/user/qrh/03-six-sevenths/routes/r4")
import model as M
_orig = np.linspace
def run(xmax, over={}):
    def lin(a, b, n, *args, **kw):
        if a == 0 and b == 0.5 and n == M.DEFAULT["xgrid"]:
            return _orig(0, xmax, n)
        return _orig(a, b, n, *args, **kw)
    M.np.linspace = lin
    try:
        return M.theta(over)[0]
    finally:
        M.np.linspace = _orig
for xm in [float(a) for a in sys.argv[1:]]:
    print(xm, run(xm), flush=True)
