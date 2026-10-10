#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/edge-theta
for N0 in 1000 2000 4000; do
  nice -n 15 python3 dual_mobius.py "{\"N0\":$N0,\"Fmax\":256,\"seeds\":8,\"out\":\"dual_unmarked_N$N0.json\"}" > dual_unmarked_N$N0.out 2>&1
done
nice -n 15 python3 dual_mobius.py '{"N0":500,"L0":13,"Fmax":128,"seeds":8,"out":"dual_marked_N500_L13.json"}' > dual_marked_N500_L13.out 2>&1
nice -n 15 python3 dual_mobius.py '{"N0":500,"L0":31,"Fmax":64,"seeds":8,"out":"dual_marked_N500_L31.json"}' > dual_marked_N500_L31.out 2>&1
nice -n 15 python3 dual_mobius.py '{"N0":1000,"L0":31,"Fmax":32,"seeds":8,"out":"dual_marked_N1000_L31.json"}' > dual_marked_N1000_L31.out 2>&1
