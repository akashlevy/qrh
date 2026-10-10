import sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/private/tmp/claude-501/-Users-akashlevy-Documents-riemann/eab6979c-f095-4107-a475-7e51d88a9fe7/scratchpad/p2")
import model_r4 as M, tsm_dropq as T
base = dict(M.DEFAULT); base.update(dgrid=321)
for name, mode, over in [("LSmarked, P^2/Y kept (D1 only, no LoA6)", "LSmarked", {}),
                         ("no split (paper Gram), P^2/Y kept", "none", {}),
                         ("no split, P^2/Y dropped", "none", {"dropq": True}),
                         ("LSmarked, P^2/Y dropped, floor 51/100, perfect counts R=1-delta", "LSmarked", {"dropq": True, "perfect": True})]:
    cfg = dict(base); cfg.update(over)
    sg, g = T.lp_for_S(0.0, cfg, mode, return_geom=True)
    print("S=0 | %-62s theta=%.6f  lx=%.4f ly=%.4f" % (name, sg, g["lx"], g["ly"]) if g else (name, sg))
