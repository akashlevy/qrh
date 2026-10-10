# matched test function: Vtilde(s) = G(1/2+s, l=1) = Gamma_C(s+1/4)Gamma_C(s+1/2)Gamma_C(s+3/4)
#   => V(x) = 8 (2pi)^{-3/2} * G^{3,0}_{0,3}((2pi)^3 x | 1/4,1/2,3/4)      (Gamma_C(s)=2(2pi)^{-s}Gamma(s))
import mpmath as mp, numpy as np
mp.mp.dps=30
def V_exact(x):
    z=(2*mp.pi)**3*x
    return 8*(2*mp.pi)**(-1.5)*mp.meijerg([[],[]],[[0.25,0.5,0.75],[]],z)
# check Mellin transform numerically at s=1: int V(x) x^{s-1} dx = Gamma_C(5/4)Gamma_C(3/2)Gamma_C(7/4)
def GC(s): return 2*(2*mp.pi)**(-s)*mp.gamma(s)
if __name__=="__main__":
    s=mp.mpf(1)
    I=mp.quad(lambda x: V_exact(x)*x**(s-1),[0,0.001,0.01,0.1,1,10])
    print("Mellin check:",I, GC(s+0.25)*GC(s+0.5)*GC(s+0.75))
    lx=np.linspace(-14,4,721)
    vals=[float(mp.log(V_exact(mp.e**t))) for t in lx]
    np.save('Vmatched_logtable.npy',np.array([lx,vals]))
    for t,v in list(zip(lx,vals))[::40]: print("log x=%6.2f  log V=%10.4f"%(t,v))
