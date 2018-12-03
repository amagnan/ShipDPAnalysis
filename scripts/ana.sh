#!/bin/sh
source /afs/cern.ch/user/t/takmete/.bashrc
source /afs/cern.ch/user/t/takmete/SHiPBuild/FairShipRun/config.sh
set -ux
echo "Starting script."
if  ! [ -e "$DP"/ana/dat/"$1"_Ana_rate1.dat ]; then
python /afs/cern.ch/user/t/takmete/SHiPBuild/FairShip/macro/Ana.py  -f /eos/experiment/ship/data/DarkPhoton/PBC-June-3/rec/"$1"_rec.root -g /eos/experiment/ship/data/DarkPhoton/PBC-June-3/geofile_full.conical.Pythia8-TGeant4.root
xrdcp "$1"_ana.root root://eospublic.cern.ch/"$DP"/ana/180112/"$1"_ana.root  
 
xrdcp "$1"_Ana_all.dat  root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_all.dat  
xrdcp "$1"_Ana_e.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_e.dat    
xrdcp "$1"_Ana_mu.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_mu.dat   
xrdcp "$1"_Ana_pi.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_pi.dat   
xrdcp "$1"_Ana_ka.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_ka.dat   
xrdcp "$1"_Ana_sum.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_sum.dat  
xrdcp "$1"_Ana_rate1.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_rate1.dat
xrdcp "$1"_Ana_rate2.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_rate2.dat
xrdcp "$1"_Ana_tau.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_tau.dat  
xrdcp "$1"_Ana_mix.dat root://eospublic.cern.ch/"$DP"/ana/180112/dat/"$1"_Ana_mix.dat  
fi
