import sys; sys.path.insert(0,'.')
from kub import *
sig=0.5; V34=math.sqrt(2*math.pi)*sig*math.exp(sig*sig*0.5625/2)
Xs=10**np.arange(4.0,5.21,0.3)
def tau_est(Nn,Cn,Kn,cl):
    m=Kn==cl; S=smooth_sums(Nn[m],Cn[m],Xs,sig)/(V34*Xs**0.75); return S
# base residues for nu = lambda^-2 (j=-2,e=0)
N0,C0,K0=coeffs_unit(-2,0,0)
t0={cl:tau_est(N0,C0,K0,cl)[-2:].mean() for cl in (True,False)}
print("tau(lambda^-2; 1) ~",t0[True]," tau(lambda^-2; 1+l3) ~",t0[False])
for pi0 in [primary((1,2)),primary((1,-2)),primary((3,2)),primary((3,-2)),primary((1,4)),primary((5,2))]:
    q=gnorm(pi0); a,b=pi0
    pi_cls1 = ((a-1)%4==0 and b%4==0)
    g_lm2 = Gj(pi0,1)*(1j)**(sym_prime((1,1),pi0)*2 % 4)   # g4(lambda^-2,pi) = (lambda^2/pi)_4 g4(1,pi)
    print("pi0",pi0,"q",q,"pi0 class1:",pi_cls1)
    for k in range(4):
        Nn,Cn,Kn=coeffs(pi0,k,j=-2,e=0,omega=0)
        row=[]
        for cl in (True,False):
            S=tau_est(Nn,Cn,Kn,cl)
            pred=None
            if k==3: pred=0
            if k==2:
                # Suzuki: tau(m^2 nu; v) = (-1)^{C(m,v)} g4(nu,m)/N(m)^{3/4} tau(nu; m v)
                Cmv = ((q-1)//4 % 2) * (0 if cl else ((5-1)//4)%2)  # N(v): v=1 -> N=1 ; v=1+lambda^3=-1+2i (N=5) -> (5-1)/4=1
                cl_mv = cl if pi_cls1 else (not cl)
                pred=(-1)**Cmv * g_lm2/q**0.75 * t0[cl_mv]
            row.append((S,pred))
        print("  k=%d:"%k, " | ".join("class%s: est %s pred %s"%('1' if i==0 else 'L', " ".join("%6.3f%+6.3fi"%(z.real,z.imag) for z in S[-3:]), ("%6.3f%+6.3fi"%(pred.real,pred.imag) if pred is not None and pred!=0 else str(pred))) for i,(S,pred) in enumerate(row)))
    sys.stdout.flush()
