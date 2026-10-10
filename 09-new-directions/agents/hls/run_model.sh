#!/bin/bash
# hls: exponent-model runs (model_hls.py = wall-spectral model_ws.py + knob hls_a)
cd /home/user/qrh/09-new-directions/agents/hls
for c in '{"hls_a":1.0,"Hoff":0.0,"bias_relax":false}' '{"hls_a":0.25,"Hoff":0.0,"bias_relax":false}' '{"hls_a":0.5,"Hoff":0.0,"bias_relax":false}' '{"hls_a":0.75,"Hoff":0.0,"bias_relax":false}' '{"hls_a":0.9,"Hoff":0.0,"bias_relax":false}' '{"hls_a":0.0,"Hoff":0.07,"bias_relax":false}' '{"hls_a":0.5,"Hoff":0.07,"bias_relax":false}' '{"hls_a":0.0,"Hoff":0.5,"bias_relax":false}' '{"hls_a":0.0,"Hoff":1.0,"bias_relax":false}'; do
  nice python3 model_hls.py "$c"
done
