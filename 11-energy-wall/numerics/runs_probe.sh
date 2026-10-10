#!/bin/bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
cd /home/user/qrh/11-energy-wall/numerics
timeout 1200 python3 probe.py '{"N0": 50000, "Fs": [1, 2, 4, 8, 16, 32], "Ymax": 262144}' > logs/probe_N50000.log 2>&1
timeout 1200 python3 probe.py '{"N0": 100000, "Fs": [2, 4, 8, 16, 32, 64], "Ymax": 262144}' > logs/probe_N100000.log 2>&1
timeout 1200 python3 probe.py '{"N0": 200000, "Fs": [4, 8, 16, 32, 64, 128], "Ymax": 262144}' > logs/probe_N200000.log 2>&1
