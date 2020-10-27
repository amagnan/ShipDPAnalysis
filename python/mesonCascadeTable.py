from array import array
import os,sys,getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:", ["date="]) 
except getopt.GetoptError:
    sys.exit()

for o,a in opts:
    if o in ('-d','--date',): date = a

prods=['meson_pi0','meson_omega','meson_eta','meson_eta1']
for prod in prods:
    inp="../data/"+date+"/"+prod+"_weight.dat"
    f=open(inp,'r')
    k=f.readlines()
    Nr, Nc =0., 0.
    for x in k:
        x=x.replace("\n","")
        x=x.split(" ")
        if float(x[7]):
            Nr+=float(x[4])
            Nc+=float(x[12])
    print prod, Nc/Nr
f.close()

