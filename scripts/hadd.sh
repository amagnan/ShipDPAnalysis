source /afs/cern.ch/user/t/takmete/.bashrc
source /afs/cern.ch/user/t/takmete/SHiPBuild/FairShipRun/config.sh
set -ux
echo "Starting script."
hadd "$1"_mass"$2"_eps"$3" "$DP"/sim/runs/"$1"_mass"$2"_eps"$3"_*.root
xrdcp "$1"_mass"$2"_eps"$3".root root://eospublic.cern.ch/"$DP"/sim/runs/"$1"_mass"$2"_eps"$3".root 
