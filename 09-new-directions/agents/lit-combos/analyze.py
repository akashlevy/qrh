"""Read runs/*.jsonl, build iso-theta threshold tables by piecewise-linear interpolation."""
import json, glob, sys
from fractions import Fraction
TARGETS = [("27/31", 27/31), ("20/23", 20/23), ("13/15", 13/15), ("6/7", 6/7)]
IGN = ("caseA", "delta0")
def load(pattern="runs/*.jsonl"):
    rows = []
    for f in sorted(glob.glob(pattern)):
        if f.endswith("test.jsonl"): continue
        for l in open(f):
            r = json.loads(l)
            if "theta" in r: rows.append(r)
    return rows
DEF = dict(wE=0.0, gram_b=1.0, kA=1.0, w1=0.0, cP=2/9, mk=0.5, perfect=False, grh=False, drop_avg=False, w2scale=0.125)
def norm(c):
    return {k: v for k, v in c.items() if k not in IGN and not (k in DEF and abs(float(v) - float(DEF[k])) < 1e-12)}
def key(c, var):
    return json.dumps({k: v for k, v in norm(c).items() if k != var}, sort_keys=True)
def curves(rows, var, default=0.0):
    fam = {}
    for r in rows:
        c = r["cfg"]
        x = norm(c).get(var, DEF.get(var, default))
        fam.setdefault(key(c, var), {})[x] = r["theta"]
    return fam
def threshold(pts, target, increasing_knob_lowers=True):
    """pts: dict x->theta; theta decreasing in x (for wE, w1, mk, cP) -> smallest x with theta<=target."""
    sgn = 1 if increasing_knob_lowers else -1
    pts = {sgn * x: t for x, t in pts.items()}
    r = _thr(pts, target)
    return (r[0], sgn * r[1])
def _thr(pts, target):
    xs = sorted(pts)
    if pts[xs[0]] <= target: return ("<=", xs[0])
    for a, b in zip(xs, xs[1:]):
        ta, tb = pts[a], pts[b]
        if tb <= target < ta:
            return ("~", a + (ta - target) * (b - a) / (ta - tb))
    return (">", xs[-1])
def fmt(t):
    s, x = t
    if s == "<=": return f"free({x:g})"
    if s == ">": return f"beyond({x:g})"
    return f"{x:.4f}"
if __name__ == "__main__":
    rows = load()
    var = sys.argv[1] if len(sys.argv) > 1 else "wE"
    fam = curves(rows, var, default={"wE": 0.0, "gram_b": 1.0, "kA": 1.0, "w1": 0.0}.get(var, 0.0))
    for k in sorted(fam):
        pts = fam[k]
        if len(pts) < 2: continue
        inc = var not in ("gram_b", "kA")
        print(k, "n=%d" % len(pts), " | ".join(f"{n}: {fmt(threshold(pts, t, inc))}" for n, t in TARGETS),
              "| theta(base)=%.6f" % pts[min(pts) if inc else max(pts)])
