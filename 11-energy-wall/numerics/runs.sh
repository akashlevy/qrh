#!/bin/bash
# production runs for steps 2-3 (each well under 20 min; 2 BLAS/OpenMP threads)
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd /home/user/qrh/11-energy-wall/numerics
COMMON='"singles": [7, 13, 19, 31, 43, 61, 97, 151, 211, 307, 409, 601, 1009], "gram_L0": 300, "slots": [30, 100, 300, 1000], "slots2": [[31, 61]], "seeds": 4'
for N0 in "$@"; do
  python3 run_energy.py "{\"N0\": $N0, $COMMON}" > logs/run_N$N0.log 2>&1
done
