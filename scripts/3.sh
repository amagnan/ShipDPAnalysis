#!/bin/sh
#echo -e "qcd "{4,5}.{0..9}" "{2.0,4.0,6.0,8.0}e-{5..10}"\n"
source /afs/cern.ch/user/t/takmete/.bashrc
source /afs/cern.ch/user/t/takmete/SHiPBuild/FairShipRun/config.sh
set -ux
echo "Starting script."
prod=$1
mass=$2
eps=$3
python /afs/cern.ch/user/t/takmete/SHiPBuild/FairShip/macro/run_simScript.py --nEvents 5000 --epsilon "$3" --mass "$2" 
xrdcp ship.conical.Pythia8-TGeant4.root root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3".root
