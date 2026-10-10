/*
 * ffenergy.c  --  exact single-mark energies over F_q[T]  (agent function-field, 11-energy-wall)
 *
 * For q prime = 1 mod 6, rows k (monic squarefree, deg kappa), columns n (monic squarefree, deg nu),
 * marks P_j (monic irreducible, deg lambda), computes EXACTLY (in Z[omega])
 *
 *   S_{k,j} = sum_n  conj chi6(Res(n,n')) * chi2(Res(k,n)) * chi6(Res(P_j,n))^{-2}
 *
 * which equals  c_nu^{-1} * sum_n beta(n) psi_k(n) chi_{P_j}(n)^{-2}  (beta = gamma_2 conj alpha;
 * see verify_gauss.py / verify_beta.py; |c_nu| = 1, the constant drops out of |S|^2).
 * Also d_{k,j} = #{n squarefree: (n, k P_j) = 1}  (the exact n-diagonal of row k).
 *
 * usage: ffenergy infile outfile [nthreads]
 *   infile: "q kappa nu nmarks lambda" then nmarks lines with lambda+1 coefficients (low->high, monic)
 *   outfile (binary int64): q kappa nu nmarks lambda nk, then per squarefree k:
 *        kidx, then per mark: A, B (S = A + B*omega), d, coprime(k,P)
 *
 * Method: per row k, a table X_k[r] = chi2(Res(k,r)) over residues r mod k, built along orbits of
 * multiplication by a unit u = T + c (one Euclid resultant per orbit); then r = n mod k is obtained
 * digitwise from n = n_hi T^kappa + n_lo  (n mod k = n_lo + (n_hi T^kappa mod k)).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

#define MAXD 24
static int q, kappa, nu, nm, lam;
static int dlog_[64], inv_[64];
static int8_t *W;          /* W[n*nm + j]: sextic exponent of column (n, mark j), -1 if zero */
static int (*marks)[MAXD + 1];
static long ctrl_seed = 0;
static char *rowmask = NULL;

static int powmod(long b, long e) { long r = 1; b %= q; if (b < 0) b += q; while (e) { if (e & 1) r = r * b % q; b = b * b % q; e >>= 1; } return (int)r; }

/* Res(a,b) = lc(a)^{deg b} prod_{a(t)=0} b(t); a,b given as coefficient arrays with degrees da, db
   (trailing zeros trimmed by caller or here). Returns value in [0,q). */
static int resultant(const int *a0, int da, const int *b0, int db) {
    int A[MAXD + 2], B[MAXD + 2], R[MAXD + 2];
    memset(A, 0, sizeof(A)); memset(B, 0, sizeof(B));
    while (da >= 0 && a0[da] == 0) da--;
    while (db >= 0 && b0[db] == 0) db--;
    if (da < 0 || db < 0) return 0;
    memcpy(A, a0, sizeof(int) * (da + 1)); memcpy(B, b0, sizeof(int) * (db + 1));
    long res = 1;
    for (;;) {
        if (db == 0) return (int)(res * powmod(B[0], da) % q);
        if (da == 0) return (int)(res * powmod(A[0], db) % q);
        if (da < db) {               /* Res(a,b) = (-1)^{da db} Res(b,a) */
            memcpy(R, A, sizeof(A)); memcpy(A, B, sizeof(A)); memcpy(B, R, sizeof(A));
            int t = da; da = db; db = t;
            if ((da * db) & 1) res = (q - res) % q;
            continue;
        }
        int ib = inv_[B[db]];
        for (int i = da; i >= db; i--) {
            int c = (int)((long)A[i] * ib % q);
            if (c) for (int j = 0; j <= db; j++) { int v = (int)((A[i - db + j] - (long)c * B[j]) % q); A[i - db + j] = v < 0 ? v + q : v; }
        }
        int dr = db - 1; while (dr >= 0 && A[dr] == 0) dr--;
        if (dr < 0) return 0;
        /* Res(a,b) = (-1)^{da db} lc(b)^{da-dr} Res(b, a mod b) */
        if ((da * db) & 1) res = (q - res) % q;
        res = res * powmod(B[db], da - dr) % q;
        memcpy(R, A, sizeof(A));            /* remainder: entries >= db are 0 */
        memcpy(A, B, sizeof(A)); da = db;
        memcpy(B, R, sizeof(A)); db = dr;
    }
}

static void idx2poly(long idx, int d, int *p) { for (int i = 0; i < d; i++) { p[i] = (int)(idx % q); idx /= q; } p[d] = 1; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: ffenergy infile outfile [threads]\n"); return 1; }
    FILE *fi = fopen(argv[1], "r");
    if (fscanf(fi, "%d %d %d %d %d", &q, &kappa, &nu, &nm, &lam) != 5) return 2;
    marks = calloc(nm, sizeof(*marks));
    for (int j = 0; j < nm; j++) for (int i = 0; i <= lam; i++) if (fscanf(fi, "%d", &marks[j][i]) != 1) return 3;
    fclose(fi);
    if (argc > 3) omp_set_num_threads(atoi(argv[3]));
    if (getenv("FFRANDOM")) ctrl_seed = atol(getenv("FFRANDOM"));
    { int g = 0; for (int c = 2; c < q; c++) { int ok = 1; for (int e = 1; e < q - 1; e++) if (powmod(c, e) == 1) { ok = 0; break; } if (ok) { g = c; break; } }
      long x = 1; for (int j = 0; j < q - 1; j++) { dlog_[x] = j; x = x * g % q; }
      for (int c = 1; c < q; c++) inv_[c] = powmod(c, q - 2); }
    long qn = 1, qk = 1; for (int i = 0; i < nu; i++) qn *= q; for (int i = 0; i < kappa; i++) qk *= q;
    long nhi = qn / qk;
    W = (int8_t *)malloc(qn * nm);
    /* columns */
    #pragma omp parallel for schedule(static)
    for (long n = 0; n < qn; n++) {
        int p[MAXD + 1], dp[MAXD + 1];
        idx2poly(n, nu, p);
        for (int i = 1; i <= nu; i++) dp[i - 1] = (int)((long)i * p[i] % q);
        int rd = nu == 0 ? 1 : resultant(p, nu, dp, nu - 1);
        for (int j = 0; j < nm; j++) {
            if (rd == 0) { W[n * nm + j] = -1; continue; }
            int rp = resultant(marks[j], lam, p, nu);
            if (rp == 0) { W[n * nm + j] = -1; continue; }
            int e = (-(dlog_[rd] % 6) - 2 * (dlog_[rp] % 6)) % 6; if (e < 0) e += 6;
            if (ctrl_seed) {   /* random-phase control: same support, i.i.d. uniform sextic phases */
                uint64_t z = (uint64_t)ctrl_seed * 0x9E3779B97F4A7C15ULL + (uint64_t)n * 0xBF58476D1CE4E5B9ULL + (uint64_t)j * 0x94D049BB133111EBULL;
                z ^= z >> 31; z *= 0xD6E8FEB86659FD93ULL; z ^= z >> 32;
                e = (int)(z % 6);
            }
            W[n * nm + j] = (int8_t)e;
        }
    }
    if (getenv("FFROWS")) {   /* comma-separated list of row indices to process */
        rowmask = (char *)calloc(qk, 1);
        char *str = strdup(getenv("FFROWS")); char *tok = strtok(str, ",");
        while (tok) { long r = atol(tok); if (r >= 0 && r < qk) rowmask[r] = 1; tok = strtok(NULL, ","); }
    }
    int64_t *out = (int64_t *)calloc(qk * (1 + 4 * nm), sizeof(int64_t));
    char *issf = (char *)calloc(qk, 1);
    #pragma omp parallel
    {
        int8_t *X = (int8_t *)malloc(qk);
        char *vis = (char *)malloc(qk);
        int64_t (*cnt)[6] = calloc(nm, sizeof(*cnt));
        #pragma omp for schedule(dynamic, 4)
        for (long kidx = 0; kidx < qk; kidx++) {
            if (rowmask && !rowmask[kidx]) continue;
            int k[MAXD + 1], dk[MAXD + 1];
            idx2poly(kidx, kappa, k);
            for (int i = 1; i <= kappa; i++) dk[i - 1] = (int)((long)i * k[i] % q);
            if (kappa > 0 && resultant(k, kappa, dk, kappa - 1) == 0) continue;
            issf[kidx] = 1;
            /* unit u = T + c with k(-c) != 0 */
            int c = 0, ku = 0;
            for (c = 0; c < q; c++) { int u[2] = {c, 1}; ku = resultant(k, kappa, u, 1); if (ku) break; }
            int eu = 3 * (dlog_[ku] & 1);
            memset(vis, 0, qk);
            for (long r = 0; r < qk; r++) {
                if (vis[r]) continue;
                int rp[MAXD + 1]; long t = r; for (int i = 0; i < kappa; i++) { rp[i] = (int)(t % q); t /= q; }
                int R = resultant(k, kappa, rp, kappa - 1);
                if (R == 0) { X[r] = -1; vis[r] = 1; continue; }
                int e = 3 * (dlog_[R] & 1);
                long cur = r; int cd[MAXD + 1]; memcpy(cd, rp, sizeof(int) * kappa);
                do {
                    X[cur] = (int8_t)e; vis[cur] = 1;
                    /* cur <- (T + c) cur mod k */
                    int top = cd[kappa - 1], nd[MAXD + 1];
                    for (int i = 0; i < kappa; i++) {
                        long v = (i ? cd[i - 1] : 0) - (long)top * k[i] + (long)c * cd[i];
                        v %= q; if (v < 0) v += q; nd[i] = (int)v;
                    }
                    memcpy(cd, nd, sizeof(int) * kappa);
                    cur = 0; for (int i = kappa - 1; i >= 0; i--) cur = cur * q + cd[i];
                    e = (e + eu) % 6;
                } while (cur != r);
            }
            /* V_i = T^i mod k, i = kappa..nu */
            int V[MAXD + 2][MAXD + 1];
            { int cd[MAXD + 1]; for (int i = 0; i < kappa; i++) cd[i] = (int)((q - k[i]) % q); /* T^kappa mod k = -k_low */
              for (int i = kappa; i <= nu; i++) {
                  memcpy(V[i], cd, sizeof(int) * kappa);
                  int top = cd[kappa - 1], nd[MAXD + 1];
                  for (int a = 0; a < kappa; a++) { long v = (a ? cd[a - 1] : 0) - (long)top * k[a]; v %= q; if (v < 0) v += q; nd[a] = (int)v; }
                  memcpy(cd, nd, sizeof(int) * kappa);
              } }
            memset(cnt, 0, sizeof(int64_t) * 6 * nm);
            long qkq = qk / q;
            if (nu < kappa) {   /* n is its own residue: index(n mod k) = index(n) + q^nu */
                for (long n = 0; n < qn; n++) {
                    int x = X[n + qn];
                    if (x < 0) continue;
                    const int8_t *w = W + n * nm;
                    for (int j = 0; j < nm; j++) { int e = w[j]; if (e >= 0) cnt[j][(e + x) % 6]++; }
                }
            }
            for (long h = 0; h < nhi; h++) {
                int s[MAXD + 1];
                for (int i = 0; i < kappa; i++) s[i] = V[nu][i];
                long t = h;
                for (int i = kappa; i < nu; i++) { int ci = (int)(t % q); t /= q; if (ci) for (int a = 0; a < kappa; a++) s[a] = (int)((s[a] + (long)ci * V[i][a]) % q); }
                /* odometer over n_lo digits 1..kappa-1 with offsets */
                long off[MAXD + 1][64];
                long pw = 1;
                for (int i = 0; i < kappa; i++) { for (int d = 0; d < q; d++) off[i][d] = ((d + s[i]) % q) * pw; pw *= q; }
                int dig[MAXD + 1]; memset(dig, 0, sizeof(dig));
                long hi = 0; for (int i = 1; i < kappa; i++) hi += off[i][0];
                int s0 = s[0];
                for (long blk = 0; blk < qkq; blk++) {
                    const int8_t *Wn = W + (h * qk + blk * q) * nm;
                    const int8_t *Xr = X + hi;
                    for (int d0 = 0; d0 < q; d0++) {
                        int x = Xr[d0 + s0 < q ? d0 + s0 : d0 + s0 - q];
                        if (x < 0) continue;
                        const int8_t *w = Wn + (long)d0 * nm;
                        for (int j = 0; j < nm; j++) { int e = w[j]; if (e >= 0) cnt[j][(e + x) % 6]++; }
                    }
                    for (int i = 1; i < kappa; i++) {
                        hi -= off[i][dig[i]]; dig[i]++;
                        if (dig[i] < q) { hi += off[i][dig[i]]; break; }
                        dig[i] = 0; hi += off[i][0];
                    }
                }
            }
            int64_t *o = out + kidx * (1 + 4 * nm);
            o[0] = kidx;
            for (int j = 0; j < nm; j++) {
                int64_t *cc = cnt[j];
                o[1 + 4 * j] = cc[0] + cc[1] - cc[3] - cc[4];
                o[2 + 4 * j] = cc[1] + cc[2] - cc[4] - cc[5];
                o[3 + 4 * j] = cc[0] + cc[1] + cc[2] + cc[3] + cc[4] + cc[5];
                o[4 + 4 * j] = resultant(k, kappa, marks[j], lam) != 0;
            }
        }
        free(X); free(vis); free(cnt);
    }
    FILE *fo = fopen(argv[2], "wb");
    long nk = 0; for (long i = 0; i < qk; i++) nk += issf[i];
    int64_t hdr[6] = {q, kappa, nu, nm, lam, nk};
    fwrite(hdr, sizeof(int64_t), 6, fo);
    for (long i = 0; i < qk; i++) if (issf[i]) fwrite(out + i * (1 + 4 * nm), sizeof(int64_t), 1 + 4 * nm, fo);
    fclose(fo);
    return 0;
}
