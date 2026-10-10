#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/alt-detect
H=1.1428571428571428
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H}" > runs/c2_nopole_h87.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"caseA\": \"paper\"}" > runs/c2_nopole_h87_caseA.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": 1.375, \"cB\": 3}" > runs/c3_nopole.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": 1.5, \"cB\": 4}" > runs/c4_nopole.json &
wait
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"grh\": true}" > runs/c2_nopole_grh.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"perfect\": true}" > runs/c2_nopole_perfect.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"cP\": 0.0, \"mk\": 0.0}" > runs/c2_nopole_noIP.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"gram_q\": 1.0}" > runs/c2_nopole_gq.json &
wait
