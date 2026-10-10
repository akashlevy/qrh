import numpy as np, sys
sys.path.insert(0,"/Users/akashlevy/Documents/riemann/gram_split/wf2")
from a2lib import *
np.set_printoptions(linewidth=220, precision=1, suppress=True)
for Q in (3000,6000):
    R=Rows(Q); R.set_scales(list(range(12)))
    sf=np.where(R.sf)[0]; Om=R.Om[sf]
    pi=R.V[:,sf]*R.zetaZ[sf][None,:]
    noise_sd=np.sqrt(np.array([ (Om**2*np.abs(R.V[i,sf])**2).sum() for i in range(12)]))
    print("Q",Q,"per-scale sd scale (sqrt sum Om^2|V_i|^2):",noise_sd)
    print("  rho f=1 sum per scale:", (R.V[:,sf]*Om).sum(1).real)
    s=(pi*Om).sum(1); print("  pi  f=1 sum per scale: Re",s.real," Im",s.imag)
    # conjugation symmetry
    key={(int(a),int(b)):t for t,(a,b) in enumerate(zip(R.ma,R.mb))}
    pairs=[(t,key[(int(a-b),int(-b))]) for t,(a,b) in enumerate(zip(R.ma,R.mb)) if (int(a-b),int(-b)) in key]
    t1=np.array([p[0] for p in pairs]); t2=np.array([p[1] for p in pairs])
    print("  rows with conjugate present: %d of %d; max|V(mbar)-conj V(m)| = %.2e"%(len(pairs),R.n,np.abs(R.Vall[:,t2]-np.conj(R.Vall[:,t1])).max()))
    both=R.sf[t1]&R.sf[t2]
    z1=R.zetaZ[t1[both]]; z2=R.zetaZ[t2[both]]
    rr=z2/np.conj(z1)
    print("  zeta(mbar)/conj(zeta(m)) distinct values:", np.unique(np.round(np.angle(rr)/(np.pi/6)).astype(int)%12))
