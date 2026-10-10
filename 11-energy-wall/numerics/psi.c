/* psi.c -- quadratic row symbols psi_k(n) = (n/k)_2 over Z[w], exact, for the energy-wall numerics.
 *
 * Each row k (squarefree, primary, prime to 6) is given as a list of factors (q, r):
 *   r >= 0 : a split prime pi | k of norm q with w = r mod pi;  (n/pi)_2 = Legendre(A + B r mod q, q)
 *   r == -1: a rational prime q with q | k (inert q, or both conjugates of a split q);
 *            (n/q)_2 = Legendre(N(n) mod q, q)   [inert: x^{(q^2-1)/2} = N(x)^{(q-1)/2};
 *                                                 split pair: (A+Br)(A+B(-1-r)) = N(n) mod q]
 * QR bit tables for every prime in the factor lists are built once (build_qr).
 * build_psi fills out[i*ncols + j] = psi_{k_i}(n_j) in {-1, 0, 1} as doubles.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

void build_qr(int64_t nq, const int64_t *qs, const int64_t *off, uint8_t *bits) {
    #pragma omp parallel for schedule(dynamic, 16)
    for (int64_t i = 0; i < nq; i++) {
        int64_t q = qs[i];
        uint8_t *b = bits + off[i];
        memset(b, 0, (size_t)((q >> 3) + 1));
        for (int64_t x = 1; 2 * x < q; x++) {
            int64_t s = (x * x) % q;
            b[s >> 3] |= (uint8_t)(1u << (s & 7));
        }
    }
}

static inline int64_t fmodq(int64_t x, int64_t q, double inv) {
    int64_t t = x - q * (int64_t)floor((double)x * inv);
    while (t < 0) t += q;
    while (t >= q) t -= q;
    return t;
}

/* rows: rowptr[nrows+1] into factor arrays fq, fr, fqi (index of q's QR table) */
void build_psi(int64_t ncols, const int64_t *A, const int64_t *B, const int64_t *N,
               int64_t nrows, const int64_t *rowptr, const int64_t *fq, const int64_t *fr, const int64_t *fqi,
               const uint8_t *bits, const int64_t *off, double *out, int nth) {
    omp_set_num_threads(nth);
    #pragma omp parallel for schedule(dynamic, 4)
    for (int64_t i = 0; i < nrows; i++) {
        double *o = out + i * ncols;
        for (int64_t j = 0; j < ncols; j++) o[j] = 1.0;
        for (int64_t f = rowptr[i]; f < rowptr[i + 1]; f++) {
            int64_t q = fq[f], r = fr[f];
            const uint8_t *b = bits + off[fqi[f]];
            double inv = 1.0 / (double)q;
            if (r >= 0) {
                for (int64_t j = 0; j < ncols; j++) {
                    int64_t x = fmodq(A[j] + B[j] * r, q, inv);
                    double s = (x == 0) ? 0.0 : (((b[x >> 3] >> (x & 7)) & 1) ? 1.0 : -1.0);
                    o[j] *= s;
                }
            } else {
                for (int64_t j = 0; j < ncols; j++) {
                    int64_t x = fmodq(N[j], q, inv);
                    double s = (x == 0) ? 0.0 : (((b[x >> 3] >> (x & 7)) & 1) ? 1.0 : -1.0);
                    o[j] *= s;
                }
            }
        }
    }
}
