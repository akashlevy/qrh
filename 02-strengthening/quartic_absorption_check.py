import cmath, math
from collections import defaultdict
def is_prime(n): return n>1 and all(n%i for i in range(2,int(n**0.5)+1))
def primary(a,b):
    for _ in range(4):
        if a%2==1 and b%2==0 and (a+b)%4==1: return a,b
        a,b=-b,a
eps_cls=defaultdict(set); mu_cls=defaultdict(set); n=0
for q in range(5,4000):
    if not is_prime(q) or q%4!=1: continue
    a=next(a for a in range(1,int(q**0.5)+1) if int(round((q-a*a)**0.5))**2==q-a*a); b=int(round((q-a*a)**0.5))
    for (A,B) in [(a,b),(a,-b)]:
        A,B=primary(A,B); p=complex(A,B); r=(-A*pow(B,-1,q))%q
        tab={1:1, r:1j, q-1:-1, (q-r)%q:-1j}
        chi=lambda t: tab[pow(t,(q-1)//4,q)]
        e=[cmath.exp(2j*math.pi*(-t*B/q)) for t in range(q)]
        g1=sum(chi(t)*e[t] for t in range(1,q)); gm1=sum(chi(t).conjugate()*e[t] for t in range(1,q))
        gq=sum((1 if pow(t,(q-1)//2,q)==1 else -1)*e[t] for t in range(1,q))
        eps=gq/math.sqrt(q); gam1=g1/math.sqrt(q); gamm1=gm1/math.sqrt(q); al=p/abs(p)
        eps_cls[(A%16,B%16)].add((round(eps.real,5),round(eps.imag,5)))
        # absorption:  mu(p)*gamma_{-1}(p) / ( conj(alpha(p)) * gamma_1(p) )  should be a ray-class function
        Phi=(-1)*gamm1/(al.conjugate()*gam1)
        mu_cls[(A%16,B%16)].add((round(Phi.real,5),round(Phi.imag,5))); n+=1
print("primes:",n)
print("eps_p constant on p mod 16:", all(len(v)==1 for v in eps_cls.values()), " values:", sorted({x for v in eps_cls.values() for x in v}))
print("Phi = mu*gamma_{-1}/(conj(alpha)*gamma_1) constant on p mod 16:", all(len(v)==1 for v in mu_cls.values()),
      " #classes:",len(mu_cls), " values:", sorted({x for v in mu_cls.values() for x in v}))
