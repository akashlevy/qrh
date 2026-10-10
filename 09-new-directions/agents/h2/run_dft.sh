#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/h2
nice -n 5 python3 h2_dft.py '{"N0": 20000, "L0": 50, "nmark": 12, "jwin": [-1, -2], "out": "runs/dft_N20000_L50.json"}' > runs/dft_N20000_L50.out 2>&1
nice -n 5 python3 h2_dft.py '{"N0": 20000, "L0": 120, "nmark": 10, "jwin": [-1, -2], "out": "runs/dft_N20000_L120.json"}' > runs/dft_N20000_L120.out 2>&1
nice -n 5 python3 h2_coh.py '{"N0": 40000, "L0s": [450, 600, 780, 1000, 1300, 1700, 2200, 2900, 3700, 4800], "out": "runs/coh_N40000.json"}' > runs/coh_N40000.out 2>&1
nice -n 5 python3 h2_kernel.py '{"N0": 10000, "L0": 1000, "jK": -3, "iters": 1501, "out": "runs/kernel_N10000.json"}' > runs/kernel_N10000.out 2>&1
