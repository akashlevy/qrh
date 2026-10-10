# Exponent model of the 08 weighted-numerator argument (formulas from paper/weighted_numerator_proof.tex)
import numpy as np, sys
alpha=5/6; c=2/9; zc=17/50
def R(d_,x):
    Dx=3-(1+4*c)*x; Px=(2-4*c*x)*(1-x); J=(alpha-d_)*Dx+d_*Px
    return 1-d_+(alpha-d_)*d_*Px/(2*J)
def worst(ell,b,kappa=0.75,grid=400,verbose=False):
    lx=(1-ell-b)/2; ly=(1-ell+b)/2; h=1-lx+ell; B=11/12-ell/4
    ds=np.linspace(1/50,3/4,grid); xs=np.linspace(0,0.5,101)
    D,X=np.meshgrid(ds,xs); q=D*X; a=(1+D)/2; Rv=R(D,X); d=h
    Eold=a-B+h*(zc-1/6)-a*ly-ell/2+q*ell+d*(Rv+D/2-zc)
    Ewt=a-B+h*(zc-1/6)-ly/2-ell/2+q*ell+d*((3*Rv+1)/4-zc)-q*(2*ly-d)/(12*kappa)
    E=np.minimum(Eold,Ewt)
    i=np.unravel_index(np.argmax(E),E.shape)
    if verbose: print(f"ell={ell:.4f} b={b:.4f} worst E={E[i]:.5f} at delta={D[i]:.3f} x={X[i]:.3f} Eold={Eold[i]:.5f} Ewt={Ewt[i]:.5f} R={Rv[i]:.4f}")
    return E[i]
if __name__=="__main__":
    print(worst(0.167,43/500,verbose=True))
    for ell in [0.16,0.165,0.167,0.168,0.17,0.1773]:
        bs=np.linspace(0,0.2,81); vals=[worst(ell,b) for b in bs]; j=int(np.argmin(vals))
        worst(ell,bs[j],verbose=True)
