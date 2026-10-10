# angular twist omega=+-1 (the actual dual family carries conj(alpha)): no X^{3/4} term even for nu=lambda^-2
import sys; sys.path.insert(0,'.')
from kub import *
sig=0.5; V34=math.sqrt(2*math.pi)*sig*math.exp(sig*sig*0.5625/2)
Xs=10**np.arange(3.0,5.21,0.3)
for om in (0,1,-1):
    N_,C_,K_=coeffs_unit(-2,0,om)
    for cl in (True,False):
        m=K_==cl; S=smooth_sums(N_[m],C_[m],Xs,sig)/(V34*Xs**0.75)
        print("omega=%+d nu=lambda^-2 class %s:"%(om,'1   ' if cl else '1+l3')," ".join("%7.3f%+7.3fi"%(z.real,z.imag) for z in S))
    for pi0 in [primary((1,4))]:
        for k in (1,2,3):
            Nn,Cn,Kn=coeffs(pi0,k,j=-2,e=0,omega=om)
            S=smooth_sums(Nn,Cn,Xs,sig)/(V34*Xs**0.75)
            print("   omega=%+d nu=lambda^-2*(1+4i)^%d (all classes):"%(om,k)," ".join("%7.3f%+7.3fi"%(z.real,z.imag) for z in S))
print("log10 X:",np.round(np.log10(Xs),2))
