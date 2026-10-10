#!/bin/bash
# usage: run_model.sh '<json>' outfile
cd /home/user/qrh/03-six-sevenths/routes/r4 && python3 model.py "$1" > "$2" 2>&1
