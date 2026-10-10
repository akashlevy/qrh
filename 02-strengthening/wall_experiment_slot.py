"""Full marked energy over a slot of primes (cross terms P1 != P2 included):
   B(k) = sum_{P in slot} N(P)^{-1/2} a(P) chi_P(k) * X^{-1/2} sum_{n' sf, (n',P)=1} a(n') chi_{n'}(P)^4 chi_{n'}(k) w(N(n')/X)
   which is the paper's marked completed sum (index n = P n', twisted multiplicativity) normalized at scale X*N(P).
   Wall (paper): rows K <= X / p0.  W predicts per-row energy = diagonal  sum_P |c(P)|^2 * (per-row energy of one T_P)."""
import math, numpy as np, wall_experiment as we
from eisenstein import cx, Z6
import sys
def run(X=3000, KMAX=30000, p0=13, control=False):
    we_main_state = {}
    # reuse construction from wall_experiment by monkeypatching main's internals: simpler to recompute here
    primes=[p for p in we.eisenstein_primes(2*X+10) if p[2]!=2]
    fields={pr[1]:we.ResidueField(*pr) for pr in primes}
    pnorm={pr[1]:(pr[2] if pr[0]=="split" else pr[2]**2) for pr in primes}
    g2={p:we.gauss_prime(F,2) for p,F in fields.items()}
    plist=sorted(fields,key=lambda p:pnorm[p])
    # columns
    cols=[]
    for e in we.elements_upto(2*X):
        if we.norm(e)<=X: continue
        pe=we.primary(e)
        if pe is None or pe!=e: continue
        rem=e; fac=[]; ok=True; Nn=we.norm(e)
        for p in plist:
            if pnorm[p]>Nn: break
            F=fields[p]; red=F.red(rem)
            if (red==0 if F.kind=="split" else red==(0,0)):
                pc=(p[0]-p[1],-p[1]); prod=we.mul(rem,pc)
                rem=(prod[0]//pnorm[p],prod[1]//pnorm[p]); fac.append(p); Nn//=pnorm[p]
                red=F.red(rem)
                if (red==0 if F.kind=="split" else red==(0,0)): ok=False;break
                if Nn==1: break
        if ok and Nn==1: cols.append((e,fac))
    def coef_of(e,fac):
        g=1+0j; acc=(1,0)
        for p in fac:
            k=fields[p].sextic_exp(acc) if acc!=(1,0) else 0
            g*=g2[p]*Z6[(4*k)%6]; acc=we.mul(acc,p)
        al=cx(e)/abs(cx(e))
        return g if control else g*al.conjugate()
    coef=np.array([coef_of(e,f) for e,f in cols])
    Ncol=np.array([we.norm(e) for e,_ in cols],float); wts=np.sin(math.pi*(Ncol/X-1))**2
    rows=we.elements_upto(KMAX)
    ra=np.array([r[0] for r in rows],np.int64); rb=np.array([r[1] for r in rows],np.int64)
    rN=ra*ra-ra*rb+rb*rb; o=np.argsort(rN); ra,rb,rN=ra[o],rb[o],rN[o]
    slot=[p for p in plist if p0<=pnorm[p]<2*p0 and fields[p].kind=="split"]
    need=sorted({p for _,f in cols for p in f}|set(slot),key=lambda p:pnorm[p])
    tab={p:we.sextic_table(fields[p],ra,rb) for p in need}
    B=np.zeros(len(rows),complex); diag=0.0; perP=[]
    for P in slot:
        T=np.zeros(len(rows),complex)
        for j,(e,fac) in enumerate(cols):
            if P in fac: continue
            ex=np.zeros(len(rows),np.int16); mask=np.ones(len(rows),bool)
            for p in fac: t=tab[p]; mask&=(t>=0); ex+=t
            kk=sum(fields[p].sextic_exp(P) for p in fac)
            T+=np.where(mask,coef[j]*wts[j]*Z6[(4*kk)%6]*np.exp(1j*math.pi/3*ex),0)
        T/=math.sqrt(X)
        aP=g2[P]*((cx(P)/abs(cx(P))).conjugate() if not control else 1)
        tP=tab[P]; chiP=np.where(tP>=0,np.exp(1j*math.pi/3*tP),0)
        cP=aP/math.sqrt(pnorm[P])
        B+=cP*chiP*T
        perP.append(np.abs(cP)**2*np.abs(chiP*T)**2)
    D=np.sum(perP,axis=0)              # pure diagonal-in-P energy per row
    EB=np.cumsum(np.abs(B)**2); ED=np.cumsum(D)
    wall=X/p0
    out=[]
    for f in (0.5,1,2,4,8,16,32,64,100):
        K=wall*f
        if K>KMAX: continue
        i=np.searchsorted(rN,K,side="right")
        out.append(f"{f:g}: {EB[i-1]/i:.4f} (diag {ED[i-1]/i:.4f})")
    print(f"slot norms [{p0},{2*p0}) ({len(slot)} primes), wall K={wall:.0f}:  per-row energy  " + ";  ".join(out), flush=True)
if __name__=="__main__":
    for p0 in (13,31,61):
        run(p0=p0, control="--control" in sys.argv)
