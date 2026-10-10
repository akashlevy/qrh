import sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf"); sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
src = open("/Users/akashlevy/Documents/riemann/gram_split/wf/lp_const_saving.py").read().split("for s0 in")[0]
exec(src)
for s0 in [0.043, 0.06, 0.1, 0.3, 1.0]:
    sg, g = th("const", s0)
    print("const s0=%.3f (P^2/Y dropped): theta = %.6f" % (s0, sg), {k: round(float(v), 4) for k, v in g.items()}, flush=True)
