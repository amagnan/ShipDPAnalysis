#!/bin/sh
source /afs/cern.ch/user/t/takmete/.bashrc
source /afs/cern.ch/user/t/takmete/SHiPBuild/FairShipRun/config.sh
set -ux
echo "Starting script."

python /afs/cern.ch/user/t/takmete/myAna.py  -f /eos/experiment/ship/data/DarkPhoton/PBC-June-3/rec/"$1"_rec.root -g /eos/experiment/ship/data/DarkPhoton/PBC-June-3/geofile_full.conical.Pythia8-TGeant4.root

xrdcp "$1"_ana.root root://eospublic.cern.ch/"$DP"/ana_rec/"$1"_ana.root

#xrdcp "$1"_true_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_rec_table.dat   
#xrdcp "$1"_e_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_e_table.dat     
#xrdcp "$1"_mu_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_mu_table.dat    
#xrdcp "$1"_pi_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_pi_table.dat    
#xrdcp "$1"_ka_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_ka_table.dat   
#xrdcp "$1"_2pi0_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_2pi0_table.dat
#xrdcp "$1"_3pi_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_3pi_table.dat
#xrdcp "$1"_4pi_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_4pi_table.dat
#xrdcp "$1"_tau_table.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_tau_table.dat   

xrdcp "$1"_rec.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_rec.dat
#xrdcp "$1"_e.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_e.dat  
#xrdcp "$1"_mu.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_mu.dat 
#xrdcp "$1"_pi.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_pi.dat 
#xrdcp "$1"_ka.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_ka.dat 
#xrdcp "$1"_2pi0.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_2pi0.dat
#xrdcp "$1"_3pi.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_3pi.dat
#xrdcp "$1"_4pi.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_4pi.dat
#xrdcp "$1"_tau.dat root://eospublic.cern.ch/"$DP"/ana_rec/dat/"$1"_tau.dat
