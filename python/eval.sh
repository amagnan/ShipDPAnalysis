#!/bin/bash
for mod in "$@"; do {
  python Eval_1.py -d 200922 -p "$mod" -a Rate1 &
  python Eval_1.py -d 200922 -p "$mod" -a ErrorRateM &
  python Eval_1.py -d 200922 -p "$mod" -a ErrorRateP &
} done
