#!/bin/sh
source /afs/cern.ch/user/t/takmete/.bashrc
source /afs/cern.ch/user/t/takmete/SHiPBuild/FairShipRun/config.sh
set -ux
echo "Starting script."

python /afs/cern.ch/user/t/takmete/SHiPBuild/FairShip/macro/ShipReco.py  -f /eos/experiment/ship/data/DarkPhoton/PBC-June-3/sim/"$1".root -g /eos/experiment/ship/data/DarkPhoton/PBC-June-3/geofile_full.conical.Pythia8-TGeant4.root

xrdcp "$1"_rec.root root://eospublic.cern.ch/"$DP"/rec/"$1"_rec.root
