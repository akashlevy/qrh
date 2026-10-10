import sys, numpy as np, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'.')
from rerun_model import theta, LS_dF, LS_BGL, LS_quad, make
for name,en in [("N(h)^3, BGL (AB)^{2/3}, reflected only",make(LS_BGL,3,False)),
                ("N(h)^3, de Faveri n=4, reflected only",make(LS_dF,3,False)),
                ("N(h)^3, dFDH-conjectured optimal quartic (= de Faveri n=4)",make(LS_dF,3,False)),
                ("N(h)^3, perfectly orthogonal quartic (A+B)",make(LS_quad,3,False)),
                ("N(h)^2, perfectly orthogonal (A+B)",make(LS_quad,2,False))]:
    th,p=theta(en); print(f"{name:60s} theta={th:.5f} (lx,ly,ell)={np.round(p,4)}")
