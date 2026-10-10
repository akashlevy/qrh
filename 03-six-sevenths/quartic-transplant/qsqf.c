// read primes file (a b N re im), output binary records for all squarefree primary c with N(c)<=NMAX:
//   int64 u, int64 v, int64 N, double re, double im   (gt = g4(1,c)/sqrt(Nc))
// recursion g(1,c'pi) = (c'/pi)_2 (-1)^{C(pi,c')} g(1,c') g(1,pi)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
typedef long long ll; typedef unsigned long long ull;
static ll pw(ll b, ll e, ll m){ ll r=1; b%=m; if(b<0)b+=m; while(e){ if(e&1) r=(__int128)r*b%m; b=(__int128)b*b%m; e>>=1;} return r; }
typedef struct { ll a,b,N,p,r; int inert; double complex gt; } Prime;
static int cmpP(const void*x,const void*y){ const Prime*a=x,*b=y; if(a->N!=b->N) return a->N<b->N?-1:1; if(a->a!=b->a) return a->a<b->a?-1:1; return a->b<b->b?-1:(a->b>b->b);}
// quadratic symbol (z/pi)_2 for z=u+vi: +1/-1, 0 if pi|z
static int q2(ll u, ll v, const Prime*q){
  ll p=q->p;
  if(!q->inert){ ll z=( (u%p+p)%p + (__int128)((v%p+p)%p)*q->r )%p; if(z==0) return 0; ll w=pw(z,(p-1)/2,p); return w==1?1:-1; }
  else { ll n=((__int128)u*u+(__int128)v*v)%p; if(n==0) return 0; ll w=pw(n,(p-1)/2,p); return w==1?1:-1; }
}
int main(int argc,char**argv){
  ll NMAX=atoll(argv[1]); FILE*fi=fopen(argv[2],"r"); FILE*fo=fopen(argv[3],"wb");
  int cap=1<<20; Prime*P=malloc(sizeof(Prime)*cap); int np=0; ll a,b,N; double re,im;
  while(fscanf(fi,"%lld %lld %lld %lf %lf",&a,&b,&N,&re,&im)==5){ if(N>NMAX) continue; Prime q; q.a=a;q.b=b;q.N=N;q.gt=re+I*im;
    if(b==0){ q.inert=1; q.p=-a; q.r=0; } else { q.inert=0; q.p=N; ll Bi=pw(((b%N)+N)%N,N-2,N); q.r=(__int128)((((-a)%N)+N)%N)*Bi%N; }
    P[np++]=q; }
  qsort(P,np,sizeof(Prime),cmpP); fprintf(stderr,"%d primes\n",np);
  typedef struct { ll u,v,N; int last; double complex gt; } St;
  St*stk=malloc(sizeof(St)*(np+64)*8); int sp=0; ll cnt=0;
  // write c=1
  { ll rec[3]={1,0,1}; double d[2]={1.0,0.0}; fwrite(rec,8,3,fo); fwrite(d,8,2,fo); cnt++; }
  for(int k=np-1;k>=0;k--){ St s={P[k].a,P[k].b,P[k].N,k,P[k].gt}; stk[sp++]=s; }
  while(sp>0){
    St s=stk[--sp];
    ll rec[3]={s.u,s.v,s.N}; double d[2]={creal(s.gt),cimag(s.gt)}; fwrite(rec,8,3,fo); fwrite(d,8,2,fo); cnt++;
    for(int k=s.last+1;k<np;k++){
      if(P[k].N > NMAX/s.N) break;
      Prime*q=&P[k];
      int s2=q2(s.u,s.v,q); if(s2==0) continue;
      ll Ce=(((q->N-1)/4)%2)*(((s.N-1)/4)%2);
      double sg = s2 * ((Ce%2)? -1.0: 1.0);
      ll nu=s.u*q->a - s.v*q->b, nv=s.u*q->b + s.v*q->a;
      St t={nu,nv,s.N*q->N,k,sg*s.gt*q->gt}; stk[sp++]=t;
      if(sp>(np+64)*8-2){ fprintf(stderr,"stack overflow\n"); exit(1);} 
    }
  }
  fclose(fo); fprintf(stderr,"%lld squarefree c\n",cnt); return 0;
}
