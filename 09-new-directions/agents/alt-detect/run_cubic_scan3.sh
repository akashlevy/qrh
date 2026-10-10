#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/alt-detect
H=1.1428571428571428
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"cP\": 0.0, \"mk\": 1e-9}" > runs/s3_noprimegain.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": 1.25}" > runs/s3_h125.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": 1.5}" > runs/s3_h150.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": 1.8}" > runs/s3_h180.json &
wait
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"tmax\": 1.0}" > runs/s3_tmax1.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"cP\": 0.0, \"mk\": 1e-9, \"tmax\": 1.0}" > runs/s3_worst.json &
python3 model_cubic.py "{\"pole\": false, \"hmin\": $H, \"margin_dummy\": 0, \"delta0\": 0.1}" > runs/s3_delta0.json &
python3 model_cubic.py "{\"pole\": true, \"beta_p\": 0.875, \"hmin\": 1.5}" > runs/s3_pole_h150.json &
wait
