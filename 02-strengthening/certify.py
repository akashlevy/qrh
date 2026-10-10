"""
Run the exact conditional certificate with given parameters.
Usage: python certify.py <w> <lx> <ly> <ell> <sigma0> [<D0 = detector floor delta0, default 1/50>]     (each as a fraction string, e.g. 1/5 447/1000 ...)
Prints a JSON line with the verdict.  Stdlib only.
"""
import re, sys, json, time, os
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "conditional_certificate.py")).read()
names = ["w", "lx", "ly", "ell", "sigma0", "D0"]
vals = sys.argv[1:7]
for k, v in zip(names, vals):
    num, _, den = v.partition("/")
    src = re.sub(rf"(?m)^{k}\s*=.*", f"{k} = F({num},{den or 1})", src)
src = src.replace('if __name__ == "__main__":', 'def _main():')
g = {"__name__": "certify_run"}
exec(compile(src, "conditional_certificate(params)", "exec"), g)
t0 = time.time()
ok_high = g["certify"](g["D0"], g["DMAX"], g["F"](0), g["F"](1, 2))
checks = {n: o for n, o, _ in g["checks"]}
out = dict(params=dict(zip(names, vals)), low=float(g["low"]), low_ok=g["low"] <= g["sigma0"],
           other_checks_ok=all(checks.values()), high_certified=ok_high,
           boxes=g["stats"]["boxes"], worst_box_upper=(float(g["stats"]["worst"]) if g["stats"]["worst"] is not None else None),
           seconds=round(time.time() - t0))
out["ALL_PASS"] = bool(out["low_ok"] and out["other_checks_ok"] and ok_high)
print(json.dumps(out), flush=True)
