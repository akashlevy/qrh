#!/bin/bash
# hls: second Gram scan, needed shape N < KL with small L (rows ~ 1.25 cols), growing N0 at fixed L0
cd /home/user/qrh/09-new-directions/agents/hls
while pgrep -f block_decor.py > /dev/null; do sleep 10; done
nice python3 gram_scan.py '[{"N0":10000,"K":1250,"L0":40,"quadref":false},{"N0":20000,"K":2500,"L0":40,"quadref":false},{"N0":40000,"K":5000,"L0":40,"quadref":false},{"N0":80000,"K":10000,"L0":40,"quadref":false}]' runs/scan2.json
