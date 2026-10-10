/* allprimary.c -- all primary n = A + B w (A=1 mod 3, B=0 mod 3) with N(n) <= Xmax, with tau(n) = #ideal divisors.
 * Reads <prefix>_primes.bin from cubic_data.c (needs primes up to Xmax).
 * Output <prefix>_all.bin records {int32 A, B, N, tau}.   usage: allprimary prefix Xmax */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef struct { int64_t a, b, N; int32_t kind, p; int64_t r; double gre, gim; } prime_t;
typedef struct { int32_t A, B, N, tau; } all_t;
static prime_t *PR; static int64_t NP; static all_t *OUT; static int64_t NO, CAP, XMAX;
static void dfs(int64_t start, int64_t A, int64_t B, int64_t N, int tau) {
    if (NO >= CAP) { CAP *= 2; OUT = realloc(OUT, sizeof(all_t) * CAP); }
    OUT[NO].A = (int32_t)A; OUT[NO].B = (int32_t)B; OUT[NO].N = (int32_t)N; OUT[NO].tau = tau; NO++;
    for (int64_t i = start; i < NP; i++) {
        const prime_t *P = &PR[i];
        if (N * P->N > XMAX) break;
        int64_t a = A, b = B, M = N;
        for (int k = 1; M * P->N <= XMAX; k++) {
            int64_t a2 = a * P->a - b * P->b, b2 = a * P->b + P->a * b - b * P->b;
            a = a2; b = b2; M *= P->N;
            dfs(i + 1, a, b, M, tau * (k + 1));
        }
    }
}
int main(int argc, char **argv) {
    char fn[512]; snprintf(fn, sizeof fn, "%s_primes.bin", argv[1]);
    XMAX = atoll(argv[2]);
    FILE *f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); NP = ftell(f) / sizeof(prime_t); fseek(f, 0, SEEK_SET);
    PR = malloc(sizeof(prime_t) * NP); if (fread(PR, sizeof(prime_t), NP, f) != (size_t)NP) return 1; fclose(f);
    CAP = 1 << 20; OUT = malloc(sizeof(all_t) * CAP);
    dfs(0, 1, 0, 1, 1);
    snprintf(fn, sizeof fn, "%s_all.bin", argv[1]);
    f = fopen(fn, "wb"); fwrite(OUT, sizeof(all_t), NO, f); fclose(f);
    fprintf(stderr, "# %ld primary n with N <= %ld\n", (long)NO, (long)XMAX);
    return 0;
}
