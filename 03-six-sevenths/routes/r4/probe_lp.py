import numpy as np
from energy_lp3 import emax, closed
for count in ("paper", "refined"):
    for (Mp, lp) in [(0.7303, 0.5968), (0.6806, 0.5985), (0.8, 0.4), (0.7, 0.4), (0.6, 0.3), (0.86, 0.29), (0.5, 0.5), (0.9, 0.5)]:
        e, arg = emax(Mp, lp, 0.0, count, ret_arg=True)
        print(f"{count:8s} M'={Mp} l'={lp}: LP={e:.5f} closed={closed(Mp,lp,0,count):.5f}  (1+lp-M')/2={(1+lp-Mp)/2 - 0:.4f} -> {arg}")
