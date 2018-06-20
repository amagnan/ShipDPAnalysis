#!/bin/sh

#epsilon=5e-7
#for mass in 0.3; do for eps in 1e-5 5e-6; do ./submit_ana.sh meson $mass $eps 1nh; done; done



#meson
#    for mass in 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
#pbrem
#    for mass in 0.3 0.4 0.6 0.8 1.0 1.2 1.5 2
#qcd
#    for mass in 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0

if (( "$#" != "4" ))
    then
    echo $# $*
    echo "Input parameters needed: <prod: meson, pbrem, qcd> <mass (GeV)> <epsilon> <batch queue: local, 1nh, 1nd,...>"
    exit
fi

MyProd=$1
mass=$2
eps=$3
MyQueue=$4


prod=180602
basedir=/eos/experiment/ship/data/DarkPhoton/PBC-June-3/AM/$prod/sim
localdir=/afs/cern.ch/work/a/ammagnan/DPSIM/$prod

mkdir -p $localdir

#for epsilon in 1e-4 9.5e-5 9e-5 8.5e-5 8e-5 7.5e-5 7e-5 6.5e-5 6e-5 5.5e-5 5e-5;
#for epsilon in 1e-5 9.5e-6 9e-6 8.5e-6 8e-6 7.5e-6 7e-6 6.5e-6 6e-6 5.5e-6 5e-6;
#for epsilon in 1e-6 9.5e-7 9e-7 8.5e-7 8e-7 7.5e-7 7e-7 6.5e-7 6e-7 5.5e-7 5e-7;
#for epsilon in 1e-7 9.5e-8 9e-8 8.5e-8 8e-8 7.5e-8 7e-8 6.5e-8 6e-8 5.5e-8 5e-8;
#for epsilon in 1e-8 9.5e-9 9e-9 8.5e-9 8e-9 7.5e-9 7e-9 6.5e-9 6e-9 5.5e-9 5e-9;
#for epsilon in 5e-5 4.5e-5 4e-5 3.5e-5 3e-5 2.5e-5 2e-5 1.5e-5 1e-5;
#for epsilon in 5e-6 4.5e-6 4e-6 3.5e-6 3e-6 2.5e-6 2e-6 1.5e-6 1e-6;
#for epsilon in 5e-7 4.5e-7 4e-7 3.5e-7 3e-7 2.5e-7 2e-7 1.5e-7 1e-7;
#for epsilon in 5e-8 4.5e-8 4e-8 3.5e-8 3e-8 2.5e-8 2e-8 1.5e-8 1e-8;
#for epsilon in 5e-9 4.5e-9 4e-9 3.5e-9 3e-9 2.5e-9 2e-9 1.5e-9 1e-9;
#do
epsilon=$eps
#tmpmass=`echo "$mass" | tr -d .`
#tmpeps=`echo "$eps" | tr -d 1-`
#basefile=$basedir/$MyProd/$tmpeps/$tmpmass/ship.conical.Pythia8-TGeant4.root
    basefile=$basedir/${MyProd}_mass${mass}_eps${eps}.root
    ls $basefile
    if (( "$?" == "0" )); then
	echo "file exist"
	mkdir -p $localdir/DP$mass/eps$eps/$MyProd/
	cd $localdir/DP$mass/eps$eps/$MyProd/
	echo "source ~ammagnan/shipsetup.sh" > submitAnaJob_eps${epsilon}.sh
	echo "cp ~ammagnan/workdir/myDPAna.py ." >> submitAnaJob_eps${epsilon}.sh
	echo "python myDPAna.py -f $basefile -n 0 -m $mass --epsilon $epsilon --prodmode $MyProd | tee $localdir/DP$mass/eps$eps/$MyProd/rate_${epsilon}.dat" >> submitAnaJob_eps${epsilon}.sh
    #echo "mv histAna*.root $localdir/DP$mass/eps$eps/$MyProd/" >> submitAnaJob.sh
	
	if [ "$MyQueue" != "local" ]; then
	    echo "cp -r * $localdir/DP$mass/eps$eps/$MyProd/" >> submitAnaJob_eps${epsilon}.sh
	    chmod u+x submitAnaJob_eps${epsilon}.sh
	    echo "Executing bsub -q $MyQueue submitAnaJob_eps${epsilon}.sh"
	    bsub -q $MyQueue submitAnaJob_eps${epsilon}.sh
	else
	    chmod u+x submitAnaJob_eps${epsilon}.sh
	    ./submitAnaJob_eps${epsilon}.sh
	fi
	
    fi
#done

