/* cubic_data.c -- exact normalised cubic Gauss sums over O = Z[w] and squarefree enumeration.
 *
 * Conventions (same as agents/alt-detect/conductor_test.py):
 *   primary  = a + b w with a = 1 mod 3, b = 0 mod 3  (w = exp(2 pi i/3), w^2 = -1 - w)
 *   g(1,c)   = sum_{x mod c} (x/c)_3 e(Tr(x/c)),   e(t) = exp(2 pi i t),  Tr: F -> Q
 *   G~(c)    = g(1,c) / sqrt(N c)
 *   twisted multiplicativity  g(1,ab) = (a/b)_3^2 g(1,a) g(1,b)   (a,b coprime primary)
 *
 * Phase 1: G~(pi) for every primary prime with N(pi) <= Nmax.
 *   split p = 1 mod 3: O/pi = F_p, w -> r = -a/b mod p, Tr(x/pi) = x(2a-b)/p.
 *     G(pi) = S_0 + w^e S_1 + w^{2e} S_2, S_j = sum_{k = j mod 3} cos(2 pi t g^k / p)  (real since -1 is a cube),
 *     g a primitive root, g^{(p-1)/3} = r^e.  G(conj pi) = conj G(pi).
 *     O(p) per prime: 6 independent Shoup-multiplication chains + two-level cos/sin tables.
 *   inert q = 2 mod 3 (q^2 <= Nmax): pi = -q, O/q = F_q[w], Tr(x/pi) = -(2u-v)/q.
 * Phase 2: all squarefree primary n with N(n) <= Xmax by DFS over primes (sorted by norm),
 *   G~(n pi) = G~(n) G~(pi) (n/pi)_3^2, mu.
 * Output: <out>_primes.bin  records {int64 a,b,N; int32 kind,p; int64 r; double gre,gim}
 *         <out>_sf.bin      records {int32 A,B,N,mu; double gre,gim}
 * usage: cubic_data Nmax out_prefix [nthreads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef struct { int64_t a, b, N; int32_t kind, p; int64_t r; double gre, gim; } prime_t;
typedef struct { int32_t A, B, N, mu; double gre, gim; } sf_t;

static void mulw(int64_t a1, int64_t b1, int64_t a2, int64_t b2, int64_t *a, int64_t *b) {
    *a = a1 * a2 - b1 * b2; *b = a1 * b2 + a2 * b1 - b1 * b2;
}
static int64_t md(int64_t x, int64_t m) { x %= m; return x < 0 ? x + m : x; }
static void primary(int64_t *a, int64_t *b) {
    static const int64_t U[6][2] = {{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,-1}};
    for (int i = 0; i < 6; i++) {
        int64_t x, y; mulw(*a, *b, U[i][0], U[i][1], &x, &y);
        if (md(x, 3) == 1 && md(y, 3) == 0) { *a = x; *b = y; return; }
    }
    fprintf(stderr, "no primary associate\n"); exit(1);
}
static uint64_t powmod(uint64_t x, uint64_t e, uint64_t p) {
    uint64_t r = 1; x %= p;
    while (e) { if (e & 1) r = (unsigned __int128)r * x % p; x = (unsigned __int128)x * x % p; e >>= 1; }
    return r;
}
static int64_t invmod(int64_t a, int64_t m) {
    int64_t g = m, x = 0, x1 = 1, a1 = md(a, m);
    while (a1) { int64_t q = g / a1, t = g - q * a1; g = a1; a1 = t; t = x - q * x1; x = x1; x1 = t; }
    return md(x, m);
}
static int factor_distinct(uint64_t n, uint64_t *f) {
    int k = 0;
    for (uint64_t d = 2; d * d <= n; d++) if (n % d == 0) { f[k++] = d; while (n % d == 0) n /= d; }
    if (n > 1) f[k++] = n;
    return k;
}
static uint64_t primroot(uint64_t p) {
    uint64_t f[64]; int k = factor_distinct(p - 1, f);
    for (uint64_t g = 2;; g++) {
        int ok = 1;
        for (int i = 0; i < k; i++) if (powmod(g, (p - 1) / f[i], p) == 1) { ok = 0; break; }
        if (ok) return g;
    }
}
/* F_q[w] arithmetic */
static void ffmul(int64_t u1, int64_t v1, int64_t u2, int64_t v2, int64_t q, int64_t *u, int64_t *v) {
    *u = md(u1 * u2 - v1 * v2, q); *v = md(u1 * v2 + u2 * v1 - v1 * v2, q);
}
static void ffpow(int64_t u, int64_t v, uint64_t e, int64_t q, int64_t *ru, int64_t *rv) {
    int64_t a = 1, b = 0;
    while (e) { if (e & 1) ffmul(a, b, u, v, q, &a, &b); ffmul(u, v, u, v, q, &u, &v); e >>= 1; }
    *ru = a; *rv = b;
}

/* Shoup modular multiplication by a fixed w < p < 2^31 */
static inline uint32_t shoup(uint32_t x, uint32_t w, uint32_t wp, uint32_t p) {
    uint32_t q = (uint32_t)(((uint64_t)x * wp) >> 32);
    uint32_t r = (uint32_t)(x * w - q * p);
    return r >= p ? r - p : r;
}

#define LB 11
static void gauss_split(prime_t *P) {
    const uint32_t p = (uint32_t)P->p;
    const int64_t a = P->a, b = P->b, r = P->r;
    uint64_t g = primroot(p);
    uint64_t c = powmod(g, (p - 1) / 3, p);
    int e;
    if (c == (uint64_t)r) e = 1; else if (c == (uint64_t)md(r * r, p)) e = 2;
    else { fprintf(stderr, "cube root mismatch p=%u\n", p); exit(1); }
    uint32_t t = (uint32_t)md(2 * a - b, p);
    /* two-level tables */
    uint32_t nlo = 1u << LB, nhi = (p >> LB) + 1;
    double *clo = malloc(sizeof(double) * nlo), *slo = malloc(sizeof(double) * nlo);
    double *chi = malloc(sizeof(double) * nhi), *shi = malloc(sizeof(double) * nhi);
    for (uint32_t i = 0; i < nlo; i++) { double th = 2 * M_PI * (double)i / p; clo[i] = cos(th); slo[i] = sin(th); }
    for (uint32_t i = 0; i < nhi; i++) { double th = 2 * M_PI * (double)((uint64_t)i << LB) / p; chi[i] = cos(th); shi[i] = sin(th); }
    uint32_t g6 = (uint32_t)powmod(g, 6, p);
    uint32_t g6p = (uint32_t)(((uint64_t)g6 << 32) / p);
    uint32_t y[6];
    y[0] = t;
    for (int j = 1; j < 6; j++) y[j] = (uint32_t)((uint64_t)y[j - 1] * g % p);
    double S[6] = {0, 0, 0, 0, 0, 0};
    uint32_t M = (p - 1) / 6;
    for (uint32_t k = 0; k < M; k++) {
        for (int j = 0; j < 6; j++) {
            uint32_t v = y[j], hi = v >> LB, lo = v & (nlo - 1);
            S[j] += chi[hi] * clo[lo] - shi[hi] * slo[lo];
            y[j] = shoup(v, g6, g6p, p);
        }
    }
    double S0 = S[0] + S[3], S1 = S[1] + S[4], S2 = S[2] + S[5];
    /* w^e S1 + w^{2e} S2 */
    double c1 = cos(2 * M_PI * e / 3), s1 = sin(2 * M_PI * e / 3);
    double c2 = cos(4 * M_PI * e / 3), s2 = sin(4 * M_PI * e / 3);
    double re = S0 + c1 * S1 + c2 * S2, im = s1 * S1 + s2 * S2;
    double sq = sqrt((double)p);
    P->gre = re / sq; P->gim = im / sq;
    free(clo); free(slo); free(chi); free(shi);
}

static void gauss_inert(prime_t *P) {
    int64_t q = P->p, Nq = q * q;
    /* generator of F_{q^2}^* */
    uint64_t f[64]; int k = factor_distinct((uint64_t)(Nq - 1), f);
    int64_t gu = -1, gv = -1;
    for (int64_t u = 0; u < q && gu < 0; u++) for (int64_t v = 1; v < q && gu < 0; v++) {
        int ok = 1;
        for (int i = 0; i < k; i++) { int64_t ru, rv; ffpow(u, v, (Nq - 1) / f[i], q, &ru, &rv); if (ru == 1 && rv == 0) { ok = 0; break; } }
        if (ok) { gu = u; gv = v; }
    }
    if (q == 2) { gu = 0; gv = 1; }
    int64_t cu, cv; ffpow(gu, gv, (Nq - 1) / 3, q, &cu, &cv);
    int e;
    if (cu == 0 && cv == 1) e = 1; else if (cu == q - 1 && cv == q - 1) e = 2;
    else { fprintf(stderr, "inert cube root mismatch q=%ld\n", (long)q); exit(1); }
    double *ct = malloc(sizeof(double) * q);
    for (int64_t i = 0; i < q; i++) ct[i] = cos(2 * M_PI * (double)i / q);
    double S[3] = {0, 0, 0};
    int64_t u = 1, v = 0;
    for (int64_t kk = 0; kk < Nq - 1; kk++) {
        S[kk % 3] += ct[md(-(2 * u - v), q)];
        ffmul(u, v, gu, gv, q, &u, &v);
    }
    double c1 = cos(2 * M_PI * e / 3), s1 = sin(2 * M_PI * e / 3);
    double c2 = cos(4 * M_PI * e / 3), s2 = sin(4 * M_PI * e / 3);
    P->gre = (S[0] + c1 * S[1] + c2 * S[2]) / q;
    P->gim = (s1 * S[1] + s2 * S[2]) / q;
    free(ct);
}

/* cubic symbol (n / pi)_3 as exponent j of w, n = A + B w; -1 if pi | n */
static int sym3(const prime_t *P, int64_t A, int64_t B) {
    if (P->kind == 0) {
        int64_t p = P->p, x = md(A + md(B, p) * P->r, p);
        if (x == 0) return -1;
        uint64_t v = powmod((uint64_t)x, (uint64_t)(p - 1) / 3, (uint64_t)p);
        if (v == 1) return 0;
        if (v == (uint64_t)P->r) return 1;
        return 2;
    } else {
        int64_t q = P->p, u = md(A, q), v = md(B, q), ru, rv;
        if (u == 0 && v == 0) return -1;
        ffpow(u, v, (uint64_t)(q * q - 1) / 3, q, &ru, &rv);
        if (ru == 1 && rv == 0) return 0;
        if (ru == 0 && rv == 1) return 1;
        return 2;
    }
}

static prime_t *PR; static int64_t NP;
static sf_t *SF; static int64_t NSF, CAP; static int64_t XMAX;
static const double W3re[3] = {1.0, -0.5, -0.5};
static double W3im[3];

static void dfs(int64_t start, int64_t A, int64_t B, int64_t N, int mu, double gre, double gim) {
    if (NSF >= CAP) { CAP *= 2; SF = realloc(SF, sizeof(sf_t) * CAP); }
    SF[NSF].A = (int32_t)A; SF[NSF].B = (int32_t)B; SF[NSF].N = (int32_t)N; SF[NSF].mu = mu;
    SF[NSF].gre = gre; SF[NSF].gim = gim; NSF++;
    for (int64_t i = start; i < NP; i++) {
        const prime_t *P = &PR[i];
        if (N * P->N > XMAX) break;
        int j = sym3(P, A, B);
        if (j < 0) continue;
        int j2 = (2 * j) % 3;
        /* G * G(pi) * w^{2j} */
        double xr = gre * P->gre - gim * P->gim, xi = gre * P->gim + gim * P->gre;
        double nr = xr * W3re[j2] - xi * W3im[j2], ni = xr * W3im[j2] + xi * W3re[j2];
        int64_t A2, B2; mulw(A, B, P->a, P->b, &A2, &B2);
        dfs(i + 1, A2, B2, N * P->N, -mu, nr, ni);
    }
}

static int cmpN(const void *x, const void *y) {
    const prime_t *a = x, *b = y;
    if (a->N != b->N) return a->N < b->N ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : (a->b > b->b);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s Nmax out_prefix [nthreads]\n", argv[0]); return 1; }
    int64_t Nmax = atoll(argv[1]); XMAX = Nmax;
    int nth = argc > 3 ? atoi(argv[3]) : 2;
    W3im[0] = 0; W3im[1] = sqrt(3.0) / 2; W3im[2] = -sqrt(3.0) / 2;
    /* sieve */
    char *isc = calloc(Nmax + 1, 1);
    for (int64_t i = 2; i * i <= Nmax; i++) if (!isc[i]) for (int64_t j = i * i; j <= Nmax; j += i) isc[j] = 1;
    int64_t cap = 1024; PR = malloc(sizeof(prime_t) * cap); NP = 0;
    int64_t nsplit_rat = 0;
    for (int64_t p = 2; p <= Nmax; p++) {
        if (isc[p] || p == 3) continue;
        if (NP + 2 >= cap) { cap *= 2; PR = realloc(PR, sizeof(prime_t) * cap); }
        if (p % 3 == 1) {
            /* Cornacchia p = x^2 + 3 y^2 via t^2 = -3 */
            uint64_t g = primroot(p);
            int64_t rr = (int64_t)powmod(g, (p - 1) / 3, p);           /* root of r^2 + r + 1 */
            int64_t t0 = md(2 * rr + 1, p);
            int64_t r0 = p, r1 = t0, lim = (int64_t)sqrt((double)p);
            while (r1 > lim) { int64_t tt = r0 % r1; r0 = r1; r1 = tt; }
            int64_t x = r1, y2 = p - x * x, y = -1;
            if (y2 % 3 == 0) { y2 /= 3; y = (int64_t)llround(sqrt((double)y2)); if (y * y != y2) y = -1; }
            if (y < 0) { fprintf(stderr, "Cornacchia failed p=%ld\n", (long)p); return 1; }
            int64_t a = x + y, b = 2 * y; primary(&a, &b);
            prime_t P; memset(&P, 0, sizeof P);
            P.a = a; P.b = b; P.N = p; P.kind = 0; P.p = (int32_t)p;
            P.r = md(-a * invmod(b, p), p);
            if (md(P.r * P.r + P.r + 1, p) != 0 || md(a + b * P.r, p) != 0) { fprintf(stderr, "bad r p=%ld\n", (long)p); return 1; }
            PR[NP++] = P;
            prime_t Q = P; Q.a = a - b; Q.b = -b; Q.r = md(-1 - P.r, p);   /* conjugate prime */
            PR[NP++] = Q; nsplit_rat++;
        } else if (p * p <= Nmax) {
            prime_t P; memset(&P, 0, sizeof P);
            P.a = -p; P.b = 0; P.N = p * p; P.kind = 1; P.p = (int32_t)p; P.r = -1;
            PR[NP++] = P;
        }
    }
    free(isc);
    fprintf(stderr, "# %ld primary primes (%ld split rational) with N <= %ld\n", (long)NP, (long)nsplit_rat, (long)Nmax);
    double t0 = omp_get_wtime();
    /* split: compute for the first of each conjugate pair (index even in split pairs); conj for partner */
    omp_set_num_threads(nth);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int64_t i = NP - 1; i >= 0; i--) {
        if (PR[i].kind == 1) gauss_inert(&PR[i]);
        else if (i + 1 < NP && PR[i + 1].kind == 0 && PR[i + 1].p == PR[i].p) gauss_split(&PR[i]);
    }
    for (int64_t i = 0; i + 1 < NP; i++)
        if (PR[i].kind == 0 && PR[i + 1].kind == 0 && PR[i + 1].p == PR[i].p) { PR[i + 1].gre = PR[i].gre; PR[i + 1].gim = -PR[i].gim; i++; }
    fprintf(stderr, "# Gauss sums done in %.1fs\n", omp_get_wtime() - t0);
    qsort(PR, NP, sizeof(prime_t), cmpN);
    char fn[512];
    snprintf(fn, sizeof fn, "%s_primes.bin", argv[2]);
    FILE *fo = fopen(fn, "wb"); fwrite(PR, sizeof(prime_t), NP, fo); fclose(fo);
    /* phase 2 */
    t0 = omp_get_wtime();
    CAP = 1 << 20; SF = malloc(sizeof(sf_t) * CAP); NSF = 0;
    dfs(0, 1, 0, 1, 1, 1.0, 0.0);
    fprintf(stderr, "# %ld squarefree n in %.1fs\n", (long)NSF, omp_get_wtime() - t0);
    snprintf(fn, sizeof fn, "%s_sf.bin", argv[2]);
    fo = fopen(fn, "wb"); fwrite(SF, sizeof(sf_t), NSF, fo); fclose(fo);
    return 0;
}
