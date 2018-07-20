import os,sys

mass=raw_input("mass interval (ie 5,6): " )
eps=raw_input("epsilon interval (ie 5..10 for 10^-5 to 10^-10): " )
os.system('echo -e "qcd "{'+mass+'}.{0..9}" "{2.0,4.0,6.0,8.0}e-{'+eps+'}"\n">>new.txt')
os.system('condor_submit_dag job.dag')
