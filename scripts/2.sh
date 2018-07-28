#!/bin/sh
source /afs/cern.ch/user/t/takmete/.bashrc
source /afs/cern.ch/user/t/takmete/SHiPBuild/FairShipRun/config.sh
set -ux
echo "Starting script."
prod=$1
mass=$2
eps=$3
python /afs/cern.ch/user/t/takmete/SHiPBuild/FairShip/macro/trueAna.py  -f /eos/experiment/ship/data/DarkPhoton/PBC-June-3/sim/"$1"_mass"$2"_eps"$3".root -g /eos/experiment/ship/data/DarkPhoton/PBC-June-3/geofile_full.conical.Pythia8-TGeant4.root

xrdcp "$1"_mass"$2"_eps"$3"_ana.root root://eospublic.cern.ch/"$DP"/ana/"$1"_mass"$2"_eps"$3"_ana.root
xrdcp "$1"_mass"$2"_eps"$3"_true_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_true_table.dat   
xrdcp "$1"_mass"$2"_eps"$3"_e_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_e_table.dat     
xrdcp "$1"_mass"$2"_eps"$3"_mu_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_mu_table.dat    
xrdcp "$1"_mass"$2"_eps"$3"_pi_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_pi_table.dat    
xrdcp "$1"_mass"$2"_eps"$3"_ka_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_ka_table.dat    
xrdcp "$1"_mass"$2"_eps"$3"_oth_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_oth_table.dat   
xrdcp "$1"_mass"$2"_eps"$3"_pi2_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_pi2_table.dat    
xrdcp "$1"_mass"$2"_eps"$3"_ka2_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_ka2_table.dat    
xrdcp "$1"_mass"$2"_eps"$3"_oth2_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_oth2_table.dat   
xrdcp "$1"_mass"$2"_eps"$3"_mix_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_mix_table.dat   
xrdcp "$1"_mass"$2"_eps"$3"_tau_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_tau_table.dat   
xrdcp "$1"_mass"$2"_eps"$3"_single_table.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_single_table.dat
xrdcp "$1"_mass"$2"_eps"$3"_true.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_true.dat
xrdcp "$1"_mass"$2"_eps"$3"_e.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_e.dat  
xrdcp "$1"_mass"$2"_eps"$3"_mu.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_mu.dat 
xrdcp "$1"_mass"$2"_eps"$3"_pi.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_pi.dat 
xrdcp "$1"_mass"$2"_eps"$3"_ka.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_ka.dat 
xrdcp "$1"_mass"$2"_eps"$3"_oth.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_oth.dat
xrdcp "$1"_mass"$2"_eps"$3"_pi2.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_pi2.dat 
xrdcp "$1"_mass"$2"_eps"$3"_ka2.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_ka2.dat 
xrdcp "$1"_mass"$2"_eps"$3"_oth2.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_oth2.dat
xrdcp "$1"_mass"$2"_eps"$3"_mix.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_mix.dat
xrdcp "$1"_mass"$2"_eps"$3"_tau.dat root://eospublic.cern.ch/"$DP"/ana/dat/"$1"_mass"$2"_eps"$3"_tau.dat
