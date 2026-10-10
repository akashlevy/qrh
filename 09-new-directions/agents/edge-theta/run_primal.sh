#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/edge-theta
for N0 in 10000 20000 40000 80000; do
  nice -n 15 python3 primal_energy.py "{\"N0\":$N0,\"L0\":300,\"seeds\":6,\"out\":\"primal_N$N0.json\"}" > primal_N$N0.out 2>&1
done
