#!/bin/bash
for mod in "$@"; do {
  python Eval_1.py -d 200922 -p "$mod" -a Rate1 >/afs/cern.ch/user/t/takmete/ShipDPAnalysis/data/200922/"$mod"_Rate1.txt &
  python Eval_1.py -d 200922 -p "$mod" -a ErrorRateM >/afs/cern.ch/user/t/takmete/ShipDPAnalysis/data/200922/"$mod"_ErrorRateM.txt &
  python Eval_1.py -d 200922 -p "$mod" -a ErrorRateP >/afs/cern.ch/user/t/takmete/ShipDPAnalysis/data/200922/"$mod"_ErrorRateP.txt &
} done
