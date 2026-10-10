#!/bin/bash
# usage: run_cfgs.sh outfile cfg1 cfg2 ...
out=$1; shift
for c in "$@"; do
  nice -n 10 python3 model.py "$c" >> "$out"
done
