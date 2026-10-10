/* sextic_kubota.c -- exact normalised sextic Gauss sums gamma_m(c), m = 1..5, over O = Z[w] on squarefree
 * primary odd c, and the sextic symbols (pi0 / c)_6 for a list of test primes pi0.
 * Purpose (agent hls): decide which Kubota series sum_c gamma_m(c; pi0^j) N(c)^{-s} have the theta pole,
 * i.e. which local theta coefficients tau_6(p^j) vanish; the hybrid large sieve's theta-extremal needs
 * tau_6 at k^3 P^{2 or 4}.
 *
 * Conventions = P78 / agents/edge-theta/edge_lib.py:
 *   primary: a = 1 mod 3, b = 0 mod 3; (x/pi)_6 = zeta6^k  <=>  x^{(N pi - 1)/6} = (1+w)^k mod pi;
 *   gamma_m(c) = N(c)^{-1/2} sum_{x mod c} (x/c)_6^m e(x/c),  e(z) = exp(4 pi i Im z / sqrt 3);
 *   split pi = a + b w (norm p): x = t in Z, e(t/pi) = exp(-2 pi i t b / p);
 *   inert pi = -q: x = u + v w, e(x/pi) = exp(-2 pi i v / q);
 *   CRT (no reciprocity used): gamma_m(c pi) = gamma_m(c) gamma_m(pi) zeta6^{m (k1 + k2)},
 *     k1 = exponent of (pi / c)_6 = sum_{q | c} exp (pi/q)_6,  k2 = exponent of (c / pi)_6.
 * Output records (packed): int32 N, int32 mu, double g[5][2] (m = 1..5), int8 e[T]  (e = exponent of (pi0_t / c)_6, -1 if pi0_t | c)
 * usage: sextic_kubota Xmax out.bin nthreads p0_1 p0_2 ...
 *    or: sextic_kubota Xmax out.bin -PRIMES=prev.bin p0_1 ...  (reuse the prime Gauss sums stored in prev.bin's header)   (p0 = rational primes = 1 mod 3; the primary prime with smaller a is used)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef struct { int64_t a, b, N; int32_t kind, p; int64_t r; int64_t w6p[6]; double g[5][2]; } prime_t;

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
static void ffmul(int64_t u1, int64_t v1, int64_t u2, int64_t v2, int64_t q, int64_t *u, int64_t *v) {
    *u = md(u1 * u2 - v1 * v2, q); *v = md(u1 * v2 + u2 * v1 - v1 * v2, q);
}
static void ffpow(int64_t u, int64_t v, uint64_t e, int64_t q, int64_t *ru, int64_t *rv) {
    int64_t a = 1, b = 0;
    while (e) { if (e & 1) ffmul(a, b, u, v, q, &a, &b); ffmul(u, v, u, v, q, &u, &v); e >>= 1; }
    *ru = a; *rv = b;
}
static inline uint32_t shoup(uint32_t x, uint32_t w, uint32_t wp, uint32_t p) {
    uint32_t q = (uint32_t)(((uint64_t)x * wp) >> 32);
    uint32_t r = (uint32_t)(x * w - q * p);
    return r >= p ? r - p : r;
}

/* sextic exponent of (A + B w / pi)_6, -1 if pi | . */
static int sym6(const prime_t *P, int64_t A, int64_t B) {
    if (P->kind == 0) {
        int64_t p = P->p, x = md(A + md(B, p) * P->r, p);
        if (x == 0) return -1;
        uint64_t v = powmod((uint64_t)x, (uint64_t)(p - 1) / 6, (uint64_t)p);
        for (int k = 0; k < 6; k++) if ((uint64_t)P->w6p[k] == v) return k;
        fprintf(stderr, "sym6 split failure\n"); exit(1);
    } else {
        int64_t q = P->p, u = md(A, q), v = md(B, q), ru, rv;
        if (u == 0 && v == 0) return -1;
        ffpow(u, v, (uint64_t)(q * q - 1) / 6, q, &ru, &rv);
        int64_t tu = 1, tv = 0;
        for (int k = 0; k < 6; k++) { if (tu == ru && tv == rv) return k; ffmul(tu, tv, 1, 1, q, &tu, &tv); }
        fprintf(stderr, "sym6 inert failure\n"); exit(1);
    }
}

#define LB 11
static void gauss_split(prime_t *P) {
    const uint32_t p = (uint32_t)P->p;
    uint64_t g = primroot(p);
    uint64_t h = powmod(g, (p - 1) / 6, p);
    int e0 = -1;
    for (int k = 0; k < 6; k++) if ((uint64_t)P->w6p[k] == h) e0 = k;
    if (e0 < 0) { fprintf(stderr, "6th root mismatch p=%u\n", p); exit(1); }
    /* additive: exp(-2 pi i t b / p); start y = b mod p, y <- y g */
    uint32_t t = (uint32_t)md(P->b, p);
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
    double Sr[6] = {0}, Si[6] = {0};
    uint32_t M = (p - 1) / 6;
    for (uint32_t k = 0; k < M; k++) {
        for (int j = 0; j < 6; j++) {
            uint32_t v = y[j], hi = v >> LB, lo = v & (nlo - 1);
            double c = chi[hi] * clo[lo] - shi[hi] * slo[lo];
            double s = shi[hi] * clo[lo] + chi[hi] * slo[lo];
            Sr[j] += c; Si[j] -= s;           /* exp(-i theta) */
            y[j] = shoup(v, g6, g6p, p);
        }
    }
    double sq = sqrt((double)p);
    for (int m = 1; m <= 5; m++) {
        double re = 0, im = 0;
        for (int j = 0; j < 6; j++) {
            double th = 2 * M_PI * (double)((m * e0 * j) % 6) / 6.0, c = cos(th), s = sin(th);
            re += c * Sr[j] - s * Si[j]; im += s * Sr[j] + c * Si[j];
        }
        P->g[m - 1][0] = re / sq; P->g[m - 1][1] = im / sq;
    }
    free(clo); free(slo); free(chi); free(shi);
}

static void gauss_inert(prime_t *P) {
    int64_t q = P->p, Nq = q * q;
    double Sr[6] = {0}, Si[6] = {0};
    for (int64_t u = 0; u < q; u++) for (int64_t v = 0; v < q; v++) {
        if (u == 0 && v == 0) continue;
        int k = sym6(P, u, v);
        double th = -2 * M_PI * (double)v / q;
        Sr[k] += cos(th); Si[k] += sin(th);
    }
    for (int m = 1; m <= 5; m++) {
        double re = 0, im = 0;
        for (int k = 0; k < 6; k++) {
            double th = 2 * M_PI * (double)((m * k) % 6) / 6.0, c = cos(th), s = sin(th);
            re += c * Sr[k] - s * Si[k]; im += s * Sr[k] + c * Si[k];
        }
        P->g[m - 1][0] = re / q; P->g[m - 1][1] = im / q;
    }
    (void)Nq;
}

static prime_t *PR; static int64_t NP; static int64_t XMAX; static int T;
static int8_t *E0;             /* NP x T : exponent of (pi0_t / pi)_6, -1 if equal */
static FILE *FO; static int64_t NOUT;
static double Z6r[6], Z6i[6];
static int fac[64];

#pragma pack(push, 1)
typedef struct { int32_t N, mu; double g[5][2]; } rec_t;
#pragma pack(pop)

static void dfs(int64_t start, int depth, int64_t A, int64_t B, int64_t N, int mu, double g[5][2], int8_t *e) {
    rec_t R; R.N = (int32_t)N; R.mu = mu; memcpy(R.g, g, sizeof(R.g));
    fwrite(&R, sizeof R, 1, FO); fwrite(e, 1, T, FO); NOUT++;
    int8_t e2[256]; double g2[5][2];
    for (int64_t i = start; i < NP; i++) {
        const prime_t *P = &PR[i];
        if (N * P->N > XMAX) break;
        if (P->N == 4) continue;                         /* odd only */
        int k2 = sym6(P, A, B);                          /* (c / pi)_6 */
        if (k2 < 0) continue;
        int k1 = 0;
        for (int l = 0; l < depth; l++) k1 += sym6(&PR[fac[l]], P->a, P->b);   /* (pi / c)_6 */
        for (int m = 1; m <= 5; m++) {
            double xr = g[m - 1][0] * P->g[m - 1][0] - g[m - 1][1] * P->g[m - 1][1];
            double xi = g[m - 1][0] * P->g[m - 1][1] + g[m - 1][1] * P->g[m - 1][0];
            int kk = (m * (k1 + k2)) % 6;
            g2[m - 1][0] = xr * Z6r[kk] - xi * Z6i[kk];
            g2[m - 1][1] = xr * Z6i[kk] + xi * Z6r[kk];
        }
        for (int t = 0; t < T; t++) {
            int8_t et = E0[i * T + t];
            e2[t] = (e[t] < 0 || et < 0) ? -1 : (int8_t)((e[t] + et) % 6);
        }
        int64_t A2, B2; mulw(A, B, P->a, P->b, &A2, &B2);
        fac[depth] = (int)i;
        dfs(i + 1, depth + 1, A2, B2, N * P->N, -mu, g2, e2);
    }
}

static int cmpN(const void *x, const void *y) {
    const prime_t *a = x, *b = y;
    if (a->N != b->N) return a->N < b->N ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : (a->b > b->b);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s Xmax out.bin nthreads p0...\n", argv[0]); return 1; }
    int64_t Nmax = atoll(argv[1]); XMAX = Nmax;
    int nth = 1; const char *PREV = NULL;
    if (strncmp(argv[3], "-PRIMES=", 8) == 0) PREV = argv[3] + 8; else nth = atoi(argv[3]);
    T = argc - 4;
    int64_t P0[256]; for (int t = 0; t < T; t++) P0[t] = atoll(argv[4 + t]);
    for (int k = 0; k < 6; k++) { Z6r[k] = cos(2 * M_PI * k / 6); Z6i[k] = sin(2 * M_PI * k / 6); }
    char *isc = calloc(Nmax + 1, 1);
    for (int64_t i = 2; i * i <= Nmax; i++) if (!isc[i]) for (int64_t j = i * i; j <= Nmax; j += i) isc[j] = 1;
    int64_t cap = 1024; PR = malloc(sizeof(prime_t) * cap); NP = 0;
    for (int64_t p = 2; p <= Nmax; p++) {
        if (isc[p] || p == 3) continue;
        if (NP + 2 >= cap) { cap *= 2; PR = realloc(PR, sizeof(prime_t) * cap); }
        if (p % 3 == 1) {
            uint64_t g = primroot(p);
            int64_t rr = (int64_t)powmod(g, (p - 1) / 3, p);
            int64_t t0 = md(2 * rr + 1, p);
            int64_t r0 = p, r1 = t0, lim = (int64_t)sqrt((double)p);
            while (r1 > lim) { int64_t tt = r0 % r1; r0 = r1; r1 = tt; }
            int64_t x = r1, y2 = p - x * x, y = -1;
            if (y2 % 3 == 0) { y2 /= 3; y = (int64_t)llround(sqrt((double)y2)); if (y * y != y2) y = -1; }
            if (y < 0) { fprintf(stderr, "Cornacchia failed p=%ld\n", (long)p); return 1; }
            int64_t a = x + y, b = 2 * y; primary(&a, &b);
            for (int c = 0; c < 2; c++) {
                prime_t P; memset(&P, 0, sizeof P);
                if (c == 0) { P.a = a; P.b = b; } else { P.a = a - b; P.b = -b; primary(&P.a, &P.b); }
                P.N = p; P.kind = 0; P.p = (int32_t)p;
                P.r = md(-P.a * invmod(P.b, p), p);
                if (md(P.r * P.r + P.r + 1, p) != 0 || md(P.a + P.b * P.r, p) != 0) { fprintf(stderr, "bad r\n"); return 1; }
                int64_t w6 = md(-(P.r * P.r % p), p), tt = 1;
                for (int k = 0; k < 6; k++) { P.w6p[k] = tt; tt = tt * w6 % p; }
                PR[NP++] = P;
            }
        } else if (p * p <= Nmax) {
            prime_t P; memset(&P, 0, sizeof P);
            P.a = -p; P.b = 0; P.N = p * p; P.kind = 1; P.p = (int32_t)p; P.r = -1;
            PR[NP++] = P;
        }
    }
    free(isc);
    fprintf(stderr, "# %ld primary primes with N <= %ld\n", (long)NP, (long)Nmax);
    double t0 = omp_get_wtime();
    if (!PREV) {
    omp_set_num_threads(nth);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int64_t i = NP - 1; i >= 0; i--) {
        if (PR[i].kind == 1) gauss_inert(&PR[i]); else gauss_split(&PR[i]);
    }
    fprintf(stderr, "# Gauss sums done in %.1fs\n", omp_get_wtime() - t0);
    }
    qsort(PR, NP, sizeof(prime_t), cmpN);
    if (PREV) {    /* the stored table is sorted by the same cmpN: match index by index */
        FILE *fp = fopen(PREV, "rb"); int32_t h2[2];
        if (!fp || fread(h2, sizeof h2, 1, fp) != 1) { fprintf(stderr, "cannot read %s\n", PREV); return 1; }
        fseek(fp, 16L * h2[0], SEEK_CUR);
        int64_t i = 0, j = 0;
        while (i < NP && j < h2[1]) {
            int64_t rec[3]; double gg[5][2];
            if (fread(rec, sizeof rec, 1, fp) != 1 || fread(gg, sizeof gg, 1, fp) != 1) break;
            j++;
            if (rec[0] > PR[NP - 1].N) break;
            if (PR[i].N != rec[0] || PR[i].a != rec[1] || PR[i].b != rec[2]) { fprintf(stderr, "prime table mismatch at %ld\n", (long)i); return 1; }
            memcpy(PR[i].g, gg, sizeof gg); i++;
        }
        fclose(fp);
        if (i != NP) { fprintf(stderr, "only %ld / %ld primes loaded\n", (long)i, (long)NP); return 1; }
        fprintf(stderr, "# loaded %ld prime Gauss sums from %s\n", (long)i, PREV);
    }
    /* test primes: for each rational p0, the primary prime above it with smaller (a, b) after sort */
    int idx0[256];
    for (int t = 0; t < T; t++) {
        idx0[t] = -1;
        for (int64_t i = 0; i < NP; i++) if (PR[i].kind == 0 && PR[i].p == P0[t]) { idx0[t] = (int)i; break; }
        if (idx0[t] < 0) { fprintf(stderr, "test prime %ld not found\n", (long)P0[t]); return 1; }
        fprintf(stderr, "# test prime %d: %ld = N(%ld + %ld w)\n", t, (long)P0[t], (long)PR[idx0[t]].a, (long)PR[idx0[t]].b);
    }
    E0 = malloc((size_t)NP * T);
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < NP; i++)
        for (int t = 0; t < T; t++) {
            const prime_t *Q = &PR[idx0[t]];
            int k = (i == idx0[t]) ? -1 : sym6(&PR[i], Q->a, Q->b);     /* (pi0 / pi)_6 */
            E0[i * T + t] = (int8_t)k;
        }
    /* header: T, then the T test primes (a, b), then the prime table (N, a, b, g) for validation */
    FO = fopen(argv[2], "wb");
    int32_t hdr[2] = {T, (int32_t)NP};
    fwrite(hdr, sizeof hdr, 1, FO);
    for (int t = 0; t < T; t++) { int64_t ab[2] = {PR[idx0[t]].a, PR[idx0[t]].b}; fwrite(ab, sizeof ab, 1, FO); }
    for (int64_t i = 0; i < NP; i++) {
        int64_t rec[3] = {PR[i].N, PR[i].a, PR[i].b};
        fwrite(rec, sizeof rec, 1, FO); fwrite(PR[i].g, sizeof(PR[i].g), 1, FO);
    }
    t0 = omp_get_wtime();
    double g1[5][2]; for (int m = 0; m < 5; m++) { g1[m][0] = 1; g1[m][1] = 0; }
    int8_t e1[256]; for (int t = 0; t < T; t++) e1[t] = 0;
    dfs(0, 0, 1, 0, 1, 1, g1, e1);
    fclose(FO);
    fprintf(stderr, "# %ld squarefree odd primary c in %.1fs\n", (long)NOUT, omp_get_wtime() - t0);
    return 0;
}
