import sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf"); sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
src = open("/Users/akashlevy/Documents/riemann/gram_split/wf/lp_const_saving.py").read().split("for s0 in")[0]
exec(src)
for mode in ["random", "randommarked", "LSmarked"]:
    sg, g = th(mode)
    print("mode %s (P^2/Y dropped): theta = %.6f" % (mode, sg), {k: round(v, 4) for k, v in g.items()}, flush=True)
