#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/alt-detect
for b in 0.5 0.6 0.7 0.8 0.85 0.9 1.0; do
  python3 model_cubic.py "{\"beta_p\": $b}" > runs/cubic_beta_$b.json &
done
wait
python3 model_cubic.py '{"pole": false}' > runs/cubic_nopole.json &
python3 model_cubic.py '{"pole": false, "cB": 3}' > runs/cubic_nopole_cB3.json &
python3 model_cubic.py '{"perfect": true}' > runs/cubic_perfect_b0.875.json &
python3 model_cubic.py '{"pole": false, "perfect": true}' > runs/cubic_nopole_perfect.json &
wait
