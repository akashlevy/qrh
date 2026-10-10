import time
from combo_model import *
for name, lev, bins in [('08 inputs kappa 3/4 (search bins)', PAPER08, BINS_SEARCH),
                        ('08 inputs kappa 3/4 (mid bins d-step .0025)', PAPER08, BINS_MID),
                        ('ALL PROVED, kappa dynamic (mid bins)', PROVED, BINS_MID),
                        ('Collet model (weighted only, kappa dyn) (mid bins)', replace(PAPER08, kappa='dynamic', Eold=False), BINS_MID)]:
    t = time.time(); th = theta_star(lev, bins)
    m, det = margin_T(lev, th, bins, detail=True)
    print(f'{name:52s} theta*={th:.8f} b={det["b"]:.4f} ell={det["ell"]:.6f} kappa={det["kappa"]:.6f} delta={det["delta"]:.4f} x={det["x"]} ({time.time()-t:.0f}s)', flush=True)
print('Collet b0 = 0.8748866470')
