import numpy as np, sys
sys.path.insert(0,"/Users/akashlevy/Documents/riemann/gram_split/wf2")
from a2lib import *
np.set_printoptions(linewidth=220, precision=2, suppress=True)
def classes(R, rows, f):
    if f == 1: return np.zeros(len(rows), int)
    key = (R.ma[rows] % f) * f + (R.mb[rows] % f)
    return np.unique(key, return_inverse=True)[1]
for Q in (1500,3000,6000):
    R = Rows(Q)
    sf = np.where(R.sf)[0]; Om = R.Om[sf]
    for scales in (list(range(4,12)), list(range(5,12)), list(range(6,12))):
        R.set_scales(scales); W = R.w; K=len(scales)
        X = {"rho": R.V[:, sf], "pi": R.V[:, sf]*R.zetaZ[sf][None,:]}
        out=[]
        for nm, x in X.items():
            for f in (1,3,36):
                lab = classes(R, sf, f); nc = lab.max()+1
                T = np.zeros((nc,K),complex)
                for i in range(K):
                    y = x[i]*Om
                    T[:,i] = np.bincount(lab,y.real,nc)+1j*np.bincount(lab,y.imag,nc)
                D = (x*Om[None,:]).T
                r = cross_est(T,W).sum()/cross_est(D,W).sum()
                out.append("%s f=%d: %.2f"%(nm,f,r))
        print(Q, "Y from", int(R.Ys[scales[0]]), "sig2(all rows)=%.2f"%R.sig2, " | ".join(out))
    # Y-profile of the untwisted f=1 sum of rho over sf rows, all 12 scales
    R.set_scales(list(range(12)))
    T1 = (R.V[:, sf]*Om).sum(1)
    print("   f=1 untwisted rho class sum per scale (Y=256..): ", np.round(T1.real,2))
    lab = classes(R, sf, 3)
    for c in range(lab.max()+1):
        Tc = (R.V[:, sf][:, lab==c]*Om[lab==c]).sum(1)
        print("   unit class",c," Re:",np.round(Tc.real,1)," Im:",np.round(Tc.imag,1))
