#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)

parser.add_option('-p', '--process' ,    dest='proc'             , help='process: meson, pbrem, qcd'            , default='meson')
parser.add_option(      '--massList'      ,    dest='massList'              , help='mass list', default='')
parser.add_option(      '--epsList'      ,    dest='epsList'              , help='eps list', default='')
(opt, args) = parser.parse_args()



prod=180602
basedir="/afs/cern.ch/work/a/ammagnan/DPSIM/%s/"%prod
eosdir="/eos/experiment/ship/data/DarkPhoton/PBC-June-3/AM/%s/sim"%prod

#eos mkdir -p $eosdir
#os.system('eos mkdir -p %s'%eosdir)

massvec=[str(x) for x in opt.massList.split(',')]
epsvec=[str(x) for x in opt.epsList.split(',')]

print "%s process:"%opt.proc

for mass in massvec :
    for epsilon in epsvec:
        print "Mass %s epsilon %s"%(mass,epsilon)
	#hadd /tmp/ammagnan/pbremSIM.root $basedir/DP$mass/eps$epsilon/pbrem/*/ship.conical.Pythia8-TGeant4.root
        #wrapper
        scriptFile = open('./hadd.sh', 'w')
	scriptFile.write('hadd /tmp/ammagnan/SIM%s.root %s/%s_mass%s_eps%s_run*.root\n'%(opt.proc,eosdir,opt.proc,mass,epsilon))
        scriptFile.write('if (( "$?" == "0" )); then\n')
        scriptFile.write('eos cp /tmp/ammagnan/SIM%s.root %s/%s_mass%s_eps%s.root\n'%(opt.proc,eosdir,opt.proc,mass,epsilon))
        scriptFile.write('fi\n')
	scriptFile.write('rm /tmp/ammagnan/SIM%s.root\n'%opt.proc)
        scriptFile.close()

        #submit
        os.system('chmod u+rwx hadd.sh')
        #os.system('cat hadd.sh')
        os.system('./hadd.sh')
