#!/bin/bash
# radical: "truth ceiling" runs of the P78 architecture (model copied from agents/hls/model_hls.py)
cd /home/user/qrh/09-new-directions/agents/radical
for c in '{}' \
 '{"Hoff":0.0,"bias_relax":false}' \
 '{"Hoff":0.0,"bias_relax":false,"grh":true}' \
 '{"Hoff":0.0,"bias_relax":false,"caseA":"refined"}' \
 '{"Hoff":0.0,"bias_relax":false,"caseA":"refined","grh":true}' \
 '{"Hoff":0.0,"bias_relax":false,"perfect":true}' \
 '{"grh":true}' \
 '{"wE":0.064,"bias_relax":true}' ; do
  nice python3 model_radical.py "$c"
done
