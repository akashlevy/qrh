// Worker: normalized quartic Gauss sums gt(pi) = g4(1,pi)/sqrt(N pi) for primary primes, e(z)=exp(4 pi i Re z)
// usage: qprimes NMAX K k out    (handles rational primes p with (prime index) % K == k)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
typedef long long ll; typedef unsigned long long ull;
static ll pw(ll b, ll e, ll m){ ll r=1; b%=m; if(b<0)b+=m; while(e){ if(e&1) r=(__int128)r*b%m; b=(__int128)b*b%m; e>>=1;} return r; }
static void primary(ll *a, ll *b){ for(int k=0;k<4;k++){ ll A=*a,B=*b; if((A%2!=0)&&(B%2==0)&&((((A+B)%4)+4)%4==1)) return; *a=-B; *b=A; } fprintf(stderr,"primary fail\n"); exit(1); }
typedef struct { ll x,y; } F2;
static F2 f2mul(F2 u, F2 v, ll p){ F2 w; w.x=((__int128)u.x*v.x-(__int128)u.y*v.y)%p; if(w.x<0)w.x+=p; w.y=((__int128)u.x*v.y+(__int128)u.y*v.x)%p; if(w.y<0)w.y+=p; return w; }
static F2 f2pow(F2 u, ull e, ll p){ F2 r={1,0}; while(e){ if(e&1) r=f2mul(r,u,p); u=f2mul(u,u,p); e>>=1;} return r; }
static int qsym_inert(ll x, ll y, ll p){ F2 z={((x%p)+p)%p,((y%p)+p)%p}; if(z.x==0&&z.y==0) return -1; F2 w=f2pow(z,(ull)(p*p-1)/4,p);
  if(w.x==1&&w.y==0) return 0; if(w.x==0&&w.y==1) return 1; if(w.x==p-1&&w.y==0) return 2; if(w.x==0&&w.y==p-1) return 3; fprintf(stderr,"inert fail\n"); exit(1);}
int main(int argc,char**argv){
  ll NMAX=atoll(argv[1]); int K=atoi(argv[2]), kk=atoi(argv[3]); FILE*fo=fopen(argv[4],"w");
  char*comp=calloc(NMAX+1,1); for(ll i=2;i*i<=NMAX;i++) if(!comp[i]) for(ll j=i*i;j<=NMAX;j+=i) comp[j]=1;
  double *cs=NULL,*sn=NULL; ll tab=0; ll idxp=0;
  for(ll p=3;p<=NMAX;p++){
    if(comp[p]) continue;
    idxp++; if(idxp%K!=kk) continue;
    if(p%4==1){
      ll x=0; for(ll t=2;;t++){ ll w=pw(t,(p-1)/4,p); if((__int128)w*w%p==p-1){x=w;break;} }
      ll r0=p,r1=x; while((__int128)r1*r1>p){ ll t=r0%r1; r0=r1; r1=t; }
      ll a=r1,b2=p-a*a; ll b=(ll)llround(sqrt((double)b2)); if(b*b!=b2){fprintf(stderr,"corn\n");exit(1);}
      ll A=a,B=b; primary(&A,&B);
      ll Bi=pw(((B%p)+p)%p,p-2,p); ll r=(__int128)((((-A)%p)+p)%p)*Bi%p; if((__int128)r*r%p!=p-1){fprintf(stderr,"r\n");exit(1);}
      if(tab<p){ free(cs); free(sn); cs=malloc(sizeof(double)*p); sn=malloc(sizeof(double)*p); tab=p; }
      for(ll t=0;t<p;t++){ double ang=2*M_PI*(double)t/(double)p; cs[t]=cos(ang); sn[t]=sin(ang); }
      ll twoa=(((2*A)%p)+p)%p;
      ll g=2; for(;;g++){ int ok=1; ll m=p-1; for(ll f=2;f*f<=m;f++) if(m%f==0){ if(pw(g,(p-1)/f,p)==1){ok=0;break;} while(m%f==0)m/=f; } if(ok&&m>1&&pw(g,(p-1)/m,p)==1) ok=0; if(ok) break; }
      ll w=pw(g,(p-1)/4,p); int kg = (w==1?0:(w==r?1:(w==p-1?2:3)));
      double bre[4]={0,0,0,0},bim[4]={0,0,0,0};
      ull y=twoa; ull G=g, P=p;
      for(ll j=0;j<p-1;j++){ bre[j&3]+=cs[y]; bim[j&3]+=sn[y]; y=y*G%P; }
      double complex S=0; for(int j=0;j<4;j++){ int e=(int)((j*(ll)kg)%4); double complex u=(e==0?1:(e==1?I:(e==2?-1:-I))); S+=u*(bre[j]+I*bim[j]); }
      double complex gt=S/sqrt((double)p);
      fprintf(fo,"%lld %lld %lld %.16f %.16f\n",A,B,p,creal(gt),cimag(gt));
      // conjugate prime: g(1,conj pi) = chi_{conj pi}(-1) conj(g(1,pi)), chi(-1)=(-1)^{(p-1)/4}
      double sg = (((p-1)/4)%2)? -1.0: 1.0; double complex gc = sg*conj(gt);
      ll Ac=A,Bc=-B; primary(&Ac,&Bc); // conj of primary is primary
      if(Ac!=A || Bc!=-B){fprintf(stderr,"conj primary?\n");exit(1);}
      fprintf(fo,"%lld %lld %lld %.16f %.16f\n",Ac,Bc,p,creal(gc),cimag(gc));
    } else if(p%4==3 && p*p<=NMAX){
      double complex S=0;
      for(ll x=0;x<p;x++){ double complex ex=cexp(-4*M_PI*I*(double)x/(double)p); for(ll yv=0;yv<p;yv++){ if(x==0&&yv==0) continue; int k4=qsym_inert(x,yv,p); double complex u=(k4==0?1:(k4==1?I:(k4==2?-1:-I))); S+=u*ex; } }
      double complex gt=S/(double)p;
      fprintf(fo,"%lld %lld %lld %.16f %.16f\n",-p,0LL,p*p,creal(gt),cimag(gt));
    }
  }
  fclose(fo); return 0;
}
