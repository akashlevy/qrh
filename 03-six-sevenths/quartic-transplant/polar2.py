import sys; sys.path.insert(0,'.')
from kub import *
sig=0.5; V34=math.sqrt(2*math.pi)*sig*math.exp(sig*sig*0.5625/2)
Xs=10**np.arange(3.0,5.01,0.5)
print("tau(nu;class) estimates  S/(Vhat(3/4) X^{3/4}) at log10X =",np.log10(Xs))
for j in range(-2,4):
    for e in range(2):
        N_,C_,K_=coeffs_unit(j,e,0)
        for cl in (True,False):
            m=K_==cl
            S=smooth_sums(N_[m],C_[m],Xs,sig)/(V34*Xs**0.75)
            print("nu=lambda^%d i^%d class %s:"%(j,e,'1   ' if cl else '1+l3'), " ".join("%7.3f%+7.3fi"%(z.real,z.imag) for z in S))
