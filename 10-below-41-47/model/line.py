import numpy as np
from fixM import worstE
B=41/47
for ell in [0.05,0.1,0.13,0.15,0.16,0.17,0.1773]:
    M=12*(B-5/6)+2*ell
    bs=np.linspace(0,0.3,61); v=[worstE(M,ell,b,B=B) for b in bs]; j=int(np.argmin(v))
    print(f"ell={ell:.4f} M={M:.4f} M+ell={M+ell:.4f} best b={bs[j]:.3f} worstE={v[j]:+.5f}")
