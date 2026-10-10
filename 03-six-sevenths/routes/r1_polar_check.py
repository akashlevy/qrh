import sys, math
import numpy as np
from r1_fake_diagonal import setup_columns
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/strengthening")
from eisenstein import cx
for X in (2000, 6000, 15000):
    cols, fields, pnorm, need, coef, W, mu = setup_columns(X)
    al = np.array([cx(e) / abs(cx(e)) for e, _ in cols])
    g2 = coef * al                      # coef = conj(alpha) gamma_2  ->  gamma_2
    n = len(cols)
    s_plain = np.sum(g2 * W); s_alpha = np.sum(coef * W)
    print(f"X={X:6d} cols={n:5d} sqrt(cols)={math.sqrt(n):6.1f} X^(5/6)={X**(5/6):8.1f} | "
          f"|sum gamma_2 W| = {abs(s_plain):8.1f}   |sum conj(alpha) gamma_2 W| = {abs(s_alpha):6.1f}")
