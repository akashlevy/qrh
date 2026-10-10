#!/bin/bash
# theta as a function of the (H2) saving exponent theta2 = 1 - Hoff, (H1)-on-average assumed, bias piece kept
cd /home/user/qrh/09-new-directions/agents/h2
for H in 0.0 0.02 0.04 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.15 0.2; do
  python3 model_h2.py "{\"Hoff\": $H, \"bias_relax\": false}" > runs/hoff_$H.json
done
