#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/alt-detect
python3 cubic_bigbox.py '{"pole": false, "hmin": 1.375, "cB": 3, "third": true}' 1 1.5 0.6 > runs/t4_c3_third_box.json &
python3 cubic_bigbox.py '{"pole": false, "hmin": 1.375, "cB": 3, "third": true}' 0.999 6 0.9 > runs/t4_c3_third_big.json &
python3 cubic_bigbox.py '{"pole": false, "hmin": 1.5, "cB": 4, "third": true}' 0.999 6 0.9 > runs/t4_c4_third_big.json &
python3 cubic_bigbox.py '{"pole": false, "hmin": 1.1428571428571428, "cB": 2, "third": true}' 0.999 6 0.9 > runs/t4_c2_third_big.json &
wait
