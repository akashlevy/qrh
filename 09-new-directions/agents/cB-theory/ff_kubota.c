/* Function-field test of the coprime-part conductor rule (cB-theory).
   F_7[T], sextic symbol (x/c)_6 = iota(Res(c,x)), iota(3^k) = zeta6^k.
   For monic squarefree c of degree d:  g(1,c) = Leg(disc c) * G^d * iota(Res(c,c'))   (validated in ff_validate.py;
   Leg(disc c) = (-1)^{d-omega(c)} by Pellet).
   For each test h we accumulate  N[h][d][e][s] = #{c : deg c = d, sf, (c,h)=1, arg g(1,c)/G^d = zeta^e, (h/c)_6 = zeta^s}.
   Then  S_{h,j}(d) = G^d sum_{e,s} N zeta^{e - j s}  is the degree-d coefficient of the coprime part of the sextic
   Kubota series at index h^j:  sum_{c sf, (c,h)=1} g(h^j, c) x^{deg c},  since g(h^j,c) = (h/c)^{-j} g(1,c).
   Output: one line per (h, d): the 36 counts.                                                          */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define P 7
static int LOG3[P];   /* log base 3 mod 7 */
static int INV[P];
static int md(int a){ a %= P; return a < 0 ? a + P : a; }

/* polynomial as coefficient array low->high, with explicit degree */
static int deg_of(const int *a, int n){ for (int i = n; i >= 0; i--) if (a[i]) return i; return -1; }

static int ipow(int b, int e){ int r = 1; b = md(b); while (e > 0){ if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return r; }

/* standard resultant Res(f,g) over F_7; f,g arrays of length up to 16 */
static int resultant(int *f0, int m, int *g0, int n){
    int f[16], g[16], r[16];
    memcpy(f, f0, sizeof(int) * 16); memcpy(g, g0, sizeof(int) * 16);
    int sgn = 1, acc = 1;
    while (1){
        if (m < 0 || n < 0) return 0;
        if (n == 0) return acc * sgn * ipow(g[0], m) % P;   /* Res(f, const) = const^m */
        if (m == 0) return acc * sgn * ipow(f[0], n) % P;   /* Res(const, g) = const^n */
        /* Res(f,g) = (-1)^{mn} lc(g)^{m-k} Res(g, f mod g) */
        memcpy(r, f, sizeof(int) * 16);
        int lcinv = INV[g[n]];
        for (int i = m; i >= n; i--){
            int c = r[i] * lcinv % P;
            if (c){ for (int j = 0; j <= n; j++) r[i - n + j] = md(r[i - n + j] - c * g[j]); }
        }
        int k = deg_of(r, n - 1);
        if (k < 0) return 0;
        if ((m * n) & 1) sgn = -sgn;
        acc = acc * ipow(g[n], m - k) % P;
        memcpy(f, g, sizeof(int) * 16); m = n;
        memset(g, 0, sizeof(int) * 16); memcpy(g, r, sizeof(int) * (k + 1)); n = k;
    }
}

#define NH 6
/* test h: 0: h = 1 ; 1: T ; 2: T^2+1 (irreducible) ; 3: T(T-1) ; 4: T^2 ; 5: T-1 */
static long long CNT[NH][12][6][6];

int main(int argc, char **argv){
    int D = argc > 1 ? atoi(argv[1]) : 8;
    for (int k = 0, v = 1; k < 6; k++, v = v * 3 % P) LOG3[v] = k;
    for (int a = 1; a < P; a++) for (int b = 1; b < P; b++) if (a * b % P == 1) INV[a] = b;
    memset(CNT, 0, sizeof CNT);
    /* d = 0: c = 1 */
    for (int h = 0; h < NH; h++) CNT[h][0][0][0] = 1;
    for (int d = 1; d <= D; d++){
        long long total = 1; for (int i = 0; i < d; i++) total *= P;
        int c[16], cp[16];
        for (long long idx = 0; idx < total; idx++){
            memset(c, 0, sizeof c); long long t = idx;
            for (int i = 0; i < d; i++){ c[i] = t % P; t /= P; }
            c[d] = 1;
            memset(cp, 0, sizeof cp);
            for (int i = 1; i <= d; i++) cp[i - 1] = c[i] * (i % P) % P;
            int dc = deg_of(cp, d - 1);
            int R = resultant(c, d, cp, dc);
            if (R == 0) continue;                         /* not squarefree */
            /* disc = (-1)^{d(d-1)/2} R ; Legendre symbol mod 7 */
            int disc = ((d * (d - 1) / 2) & 1) ? md(-R) : R;
            int leg = (ipow(disc, 3) == 1) ? 1 : -1;
            int e = (LOG3[R] + (leg < 0 ? 3 : 0)) % 6;
            /* evaluate c at 0, 1, and c mod (T^2+1) */
            int c0 = c[0];
            int c1 = 0; for (int i = 0; i <= d; i++) c1 += c[i]; c1 %= P;
            /* c mod T^2+1: T^2 = -1 -> T^i = (-1)^{i/2} T^{i%2} */
            int r0 = 0, r1 = 0;
            for (int i = 0; i <= d; i++){ int s = ((i / 2) & 1) ? -1 : 1; if (i & 1) r1 += s * c[i]; else r0 += s * c[i]; }
            r0 = md(r0); r1 = md(r1);
            int sgnd = (d & 1) ? P - 1 : 1;                 /* (-1)^d */
            int res_h[NH];
            res_h[0] = 1;
            res_h[1] = sgnd * c0 % P;                     /* Res(c, T)   = (-1)^d c(0) */
            res_h[2] = md(r0 * r0 + r1 * r1);             /* Res(c, T^2+1) = c(i) c(-i) */
            res_h[3] = c0 * c1 % P;                       /* Res(c, T(T-1)) = c(0) c(1) */
            res_h[4] = c0 * c0 % P;                       /* Res(c, T^2) */
            res_h[5] = sgnd * c1 % P;                     /* Res(c, T-1) */
            for (int h = 0; h < NH; h++){
                if (res_h[h] == 0) continue;
                CNT[h][d][e][LOG3[res_h[h]]]++;
            }
        }
        fprintf(stderr, "d=%d done\n", d);
    }
    for (int h = 0; h < NH; h++) for (int d = 0; d <= D; d++){
        printf("%d %d", h, d);
        for (int e = 0; e < 6; e++) for (int s = 0; s < 6; s++) printf(" %lld", CNT[h][d][e][s]);
        printf("\n");
    }
    return 0;
}
