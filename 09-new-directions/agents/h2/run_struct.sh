#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/h2
nice -n 5 python3 h2_struct.py '{"N0": 20000, "L0s": [30, 100, 300, 1000, 3000], "nmark": 60, "out": "runs/struct_N20000.json"}' > runs/struct_N20000.out 2>&1
nice -n 5 python3 h2_struct.py '{"N0": 40000, "L0s": [100, 1000, 3000], "nmark": 80, "out": "runs/struct_N40000.json"}' > runs/struct_N40000.out 2>&1
