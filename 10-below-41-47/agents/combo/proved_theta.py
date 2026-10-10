"""Best theta with PROVED inputs only (and the 08 / Collet subsets)."""
import time
from combo_model import *
cases = [
 ('08 inputs, kappa = 3/4 fixed', PAPER08),
 ('08 + kappa dynamic (Collet variable kappa)', replace(PAPER08, kappa='dynamic')),
 ('Collet PR#9 model: weighted branch only, kappa dynamic', replace(PAPER08, kappa='dynamic', Eold=False)),
 ('ALL PROVED (kappa dyn, slot powers, sextic LS short plain, E_old)', PROVED),
]
for name, lev in cases:
    t = time.time()
    r = theta_best(lev, verbose=False)
    m, det = verify(r, lev)
    print(f'{name:66s} theta={r["theta"]:.7f} ell={r["ell"]:.5f} b={r["b"]:.4f} kappa={r["kappa"]:.5f} '
          f'| fine-grid margin {m:+.2e} at delta={det["delta"]:.4f} x={det["x"]:.3f} ({time.time()-t:.0f}s)', flush=True)
print('Collet b0 = 0.874886647 (ell0 = 0.1671200785); 10499/12000 =', 10499/12000)
