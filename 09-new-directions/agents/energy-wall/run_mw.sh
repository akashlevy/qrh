#!/bin/bash
cd /home/user/qrh/09-new-directions/agents/energy-wall && python3 model_wall.py "$1" > "$2" 2>&1
