from array import array
import os,sys,getopt
from decimal import Decimal

lepto= False

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:l:", ["date=","leptophilic="]) 
except getopt.GetoptError:
    sys.exit()
for o,a in opts:
    if o in ('-l','--leptophilic',): lepto = a
    if o in ('-d','--date',): date = a

prods=['meson','pbrem','qcd']

for prod in prods:

    if not lepto: 
        outp="../data/"+date+"/"+prod+"_Rate1.dat"
        inp="../data/"+date+"/"+prod+"_Ana_rate1.dat"
    if lepto:
        outp="../data/"+date+"/"+prod+"_Rate2.dat"
        inp="../data/"+date+"/"+prod+"_Ana_rate2.dat"
    
    f=open(inp,'r')
    k=f.readlines()
    l=open(outp,'w')

    for x in k:
        x=x.replace("\n","")
        x=x.split(" ")
        if Decimal(x[2]): 
            l.write('%.5E %.9E %.9E' %(Decimal(x[0]),Decimal(x[1]),Decimal(x[2])))
            l.write('\n')
    f.close()
    l.close()
