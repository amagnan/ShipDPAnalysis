#!/bin/sh
#for run in `seq 0 4`; do for mass in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do for epsilon in 5e-5 4e-5 2e-5 1e-5 8e-6 6e-6 4e-6 2e-6 1e-6 8e-7 5e-7 1e-7 8e-8 6e-8 4e-8 2e-8 1e-8; do for prod in meson pbrem; do ./submit_job.sh 180602 $prod $mass $epsilon 0 500 1nd; done; done; done; done


if (( "$#" != "7" ))
    then
    echo $# $*
    echo "Input parameters needed: <outputdir> <prod: meson, pbrem, qcd> <mass (GeV)> <epsilon> <run> <nEvts> <batch queue: local, 1nh, 1nd,...>"
    exit
fi

MySubdir=$1
MyProd=$2
MyMass=$3
MyEps=$4
MyRun=$5
NEVTS=$6
MyQueue=$7

eosdir=/eos/experiment/ship/data/DarkPhoton/PBC-June-3/AM/$MySubdir/sim

MyDir=/afs/cern.ch/work/a/ammagnan/DPSIM/$MySubdir/DP${MyMass}/eps${MyEps}/$MyProd/${MyRun}

mkdir -p $MyDir

cd $MyDir
echo "source ~ammagnan/shipsetup.sh" > submitJob.sh

echo "cp $FAIRSHIP/macro/run_simScript.py ." >> submitJob.sh

echo "mkdir OutData" >> submitJob.sh

echo "python run_simScript.py -n $NEVTS --DarkPhoton -A $MyProd --mass $MyMass --epsilon $MyEps | tee output.log" >> submitJob.sh

if [ "$MyQueue" != "local" ]; then
    echo "eos mkdir -p $eosdir" >> submitJob.sh
    echo "eos cp \$PWD/ship.conical.Pythia8-TGeant4.root $eosdir/${MyProd}_mass${MyMass}_eps${MyEps}_run${MyRun}.root" >> submitJob.sh
    echo "if (( \"\$?\" != \"0\" )); then" >> submitJob.sh
    echo " echo ll" >> submitJob.sh
    echo " echo \"--- Problem with copy of file to EOS. Keeping locally.\"" >> submitJob.sh
    echo "else"  >> submitJob.sh
    echo " rm ship.conical.Pythia8-TGeant4.root"  >> submitJob.sh
    echo "fi"  >> submitJob.sh
    echo "cp -r * $MyDir/." >> submitJob.sh
    chmod u+x submitJob.sh
    echo "Executing bsub -q $MyQueue submitJob.sh"
    bsub -q $MyQueue submitJob.sh
    #cat submitJob.sh
else
    #cat submitJob.sh
    chmod u+x submitJob.sh
    ./submitJob.sh
fi
