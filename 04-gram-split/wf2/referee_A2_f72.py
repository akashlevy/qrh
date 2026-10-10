import numpy as np, sys
sys.path.insert(0,"/Users/akashlevy/Documents/riemann/gram_split/wf2")
from a2lib import *
rng=np.random.default_rng(3)
def classes(R, rows, f):
    if f == 1: return np.zeros(len(rows), int)
    key = (R.ma[rows] % f) * f + (R.mb[rows] % f)
    return np.unique(key, return_inverse=True)[1]
def rat(R, x, Om, lab, W):
    K=x.shape[0]; nc=lab.max()+1
    T=np.zeros((nc,K),complex)
    for i in range(K):
        y=x[i]*Om; T[:,i]=np.bincount(lab,y.real,nc)+1j*np.bincount(lab,y.imag,nc)
    D=(x*Om[None,:]).T
    num=cross_est(T,W).sum(); den=cross_est(D,W).sum()
    # jackknife
    reps=[]
    for i in range(K):
        kp=[j for j in range(K) if j!=i]
        reps.append(cross_est(T[:,kp],W[kp]).sum()/cross_est(D[:,kp],W[kp]).sum())
    reps=np.array(reps); se=np.sqrt((K-1)/K*((reps-reps.mean())**2).sum())
    return num/den, se
for Q in (3000,6000):
    R=Rows(Q, scales=list(range(5,12)))
    W=R.w
    for rowsname, rows in (("sf",np.where(R.sf)[0]),("prime",np.where(R.sf&(R.nfac==1))[0])):
        Om=R.Om[rows]
        for f in (36,72):
            lab=classes(R,rows,f)
            print("Q=%d %s f=%d rows/class=%.2f"%(Q,rowsname,f,len(rows)/(lab.max()+1)))
            for nm,x in (("rho",R.V[:,rows]),("pi",R.V[:,rows]*R.zetaZ[rows][None,:])):
                line=[]
                for k in range(-6,7):
                    h=R.alpha_m[rows]**k
                    r,s=rat(R,x*h[None,:],Om,lab,W); line.append("%.2f(%.2f)"%(r,s))
                print("  %s k=-6..6:"%nm," ".join(line))
            nul=[rat(R,R.V[:,rows]*R.zetaZ[rows][None,:]*rng.choice([-1.,1.],len(rows))[None,:],Om,lab,W)[0] for _ in range(20)]
            print("  null(pi*eps): mean %.2f sd %.2f"%(np.mean(nul),np.std(nul)))
