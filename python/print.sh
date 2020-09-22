#!/bin/bash
for mod in "$@"; do {
  echo Eval_1.py -d 200922 -p "$mod" -a Rate1
  echo Eval_1.py -d 200922 -p "$mod" -a ErrorRateM
  echo Eval_1.py -d 200922 -p "$mod" -a ErrorRateP
} done
