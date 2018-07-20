#!/bin/sh
source /afs/cern.ch/user/t/takmete/.bashrc
source /afs/cern.ch/user/t/takmete/SHiPBuild/FairShipRun/config.sh
set -ux
echo "Starting script."
prod=$1
mass=$2
eps=$3
python /afs/cern.ch/user/t/takmete/ShipDPAnalysis/python/trueAna.py  -f /eos/experiment/ship/data/DarkPhoton/PBC-June-3/sim/"$1"_mass"$2"_eps"$3".root -g /eos/experiment/ship/data/DarkPhoton/PBC-June-3/geofile_full.conical.Pythia8-TGeant4.root

xrdcp sim/"$1"_mass"$2"_eps"$3".root root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3".root
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_rec.dat     root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_true.dat   
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_e.dat       root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_e.dat     
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_mu.dat      root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_mu.dat    
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_pi.dat      root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_pi.dat    
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_ka.dat      root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_ka.dat    
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_oth.dat     root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_oth.dat   
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_mix.dat     root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_mix.dat   
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_tau.dat     root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_tau.dat   
xrdcp sim/"$1"_mass"$2"_eps"$3"_table_single.dat  root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_table_single.dat
xrdcp sim/"$1"_mass"$2"_eps"$3"_rec.dat root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_true.dat
xrdcp sim/"$1"_mass"$2"_eps"$3"_e.dat   root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_e.dat  
xrdcp sim/"$1"_mass"$2"_eps"$3"_mu.dat  root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_mu.dat 
xrdcp sim/"$1"_mass"$2"_eps"$3"_pi.dat  root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_pi.dat 
xrdcp sim/"$1"_mass"$2"_eps"$3"_ka.dat  root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_ka.dat 
xrdcp sim/"$1"_mass"$2"_eps"$3"_oth.dat root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_oth.dat
xrdcp sim/"$1"_mass"$2"_eps"$3"_mix.dat root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_mix.dat
xrdcp sim/"$1"_mass"$2"_eps"$3"_tau.dat root://eospublic.cern.ch/"$DP"/sim/"$1"_mass"$2"_eps"$3"_tau.dat
