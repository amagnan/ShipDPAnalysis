#!/bin/sh

# for prod in qcd meson pbrem; do ./getOutput.sh $prod; done

if (( "$#" != "1" ))
    then
    echo $# $*
    echo "Input parameters needed: <prod: meson, pbrem, qcd>"
    exit
fi

prod=$1

grep "Total A' production rate" /afs/cern.ch/work/a/ammagnan/DPSIM/180602/DP*/*/$prod/*.dat | awk '{print $6,$7,$8,$9}' > $prod.txt

grep "Decay vertex in vessel acceptance" /afs/cern.ch/work/a/ammagnan/DPSIM/180602/DP*/*/$prod/*.dat | awk '{print $6}' > vtx_$prod.txt

paste $prod.txt vtx_$prod.txt > $prod.dat


