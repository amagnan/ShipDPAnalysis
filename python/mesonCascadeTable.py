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
        l=open(inp.replace('weight','boosted'),'w')
        k=f.readlines()
        Nr, Nc =0., 0.
        for x in k:
                x=x.replace("\n","")
                x=x.split(" ")
                if float(x[7]):
                        l.write('%s %s %.8g %.8g' %(x[0],x[1],float(x[12])/float(x[4]),float(x[13])/float(x[7])))
                        l.write('\n')
                        Nr+=float(x[7])
                        Nc+=float(x[13])
                        if float(x[1])==1e-07 or float(x[1])==1e-06: print prod.replace('meson_',''),x[0],x[1],float(x[12])/float(x[4]),float(x[13])/float(x[7])
                        #if float(x[1])==1e-06: print float(x[4]),float(x[12]),float(x[4])/float(x[12]),float(x[7]),float(x[13]),float(x[7])/float(x[13])
        print prod, Nc/Nr
f.close()

