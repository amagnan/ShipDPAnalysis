#!/bin/bash
for mod in "$@"; do {
  python Eval_1.py -d latest-greatest -p "$mod" -a Rate1 > /Users/takmete/Documents/GitHub/ShipDPAnalysis/data/latest-greatest/"$mod".txt &
  python Eval_1.py -d latest-greatest -p "$mod" -a ErrorRateM > /Users/takmete/Documents/GitHub/ShipDPAnalysis/data/latest-greatest/"$mod"M.txt &
  python Eval_1.py -d latest-greatest -p "$mod" -a ErrorRateP > /Users/takmete/Documents/GitHub/ShipDPAnalysis/data/latest-greatest/"$mod"P.txt &
} done
