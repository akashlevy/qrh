/* probe.c -- sextic Gauss-sum probe vectors A_{m,sigma}(Y) = sum_{s in sigma} coef(s) conj chi_s(-m), exact symbols.
 * E[p * nrows + m] = sextic exponent of (-m / p)_6 in 0..5, or -1 if p | m   (p = index into the prime table)
 * s given by ptr/fp (prime-table indices of its factors), class cls[s] in [0, ncls), complex coef (cre, cim).
 * out[(cls * nrows + m)] += coef(s) * zeta6^{-e(s,m)},  e = sum_{p|s} E[p][m] mod 6.
 */
#include <stdint.h>
#include <math.h>
#include <omp.h>

void probe(int64_t nrows, const int8_t *E, int64_t ns, const int64_t *ptr, const int64_t *fp, const int64_t *cls,
           const double *cre, const double *cim, int64_t ncls, double *ore, double *oim, int nth) {
    double zr[6], zi[6];
    for (int k = 0; k < 6; k++) { zr[k] = cos(2 * M_PI * k / 6.0); zi[k] = -sin(2 * M_PI * k / 6.0); }  /* zeta6^{-k} */
    omp_set_num_threads(nth);
    #pragma omp parallel
    {
        int nt = omp_get_num_threads(), id = omp_get_thread_num();
        int64_t r0 = nrows * id / nt, r1 = nrows * (id + 1) / nt;
        for (int64_t s = 0; s < ns; s++) {
            double a = cre[s], b = cim[s];
            double *Or = ore + cls[s] * nrows, *Oi = oim + cls[s] * nrows;
            int64_t f0 = ptr[s], f1 = ptr[s + 1];
            for (int64_t m = r0; m < r1; m++) {
                int e = 0, bad = 0;
                for (int64_t f = f0; f < f1; f++) {
                    int8_t v = E[fp[f] * nrows + m];
                    if (v < 0) { bad = 1; break; }
                    e += v;
                }
                if (bad) continue;
                e %= 6;
                Or[m] += a * zr[e] - b * zi[e];
                Oi[m] += a * zi[e] + b * zr[e];
            }
        }
    }
}
