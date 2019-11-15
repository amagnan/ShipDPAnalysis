from array import array
import os,sys,getopt
import math
lepto=0
modes = ['meson','pbrem','qcd']
fracs = ['e','mu','tau','pi0','nhadron','chadron','hadron','all','other']
try:
    opts, args = getopt.getopt(sys.argv[1:], "d:")

except getopt.GetoptError:
    print 'decay to all SM particles'

for o,a in opts:
    if o in ('-d',): date = a

pathW = "../data/"+date+"/"
#f0 = open("parameter.txt","r")#dosya lazim
#l0 = f0.readlines()
for mode in modes:
    for frac in fracs:
        exec('lFix_%s_%s  = open(pathW+"%sFix_%s.dat","w")'%(mode,frac,mode,frac))
        exec('%s_%s = open(pathW+"%s_%s.dat","r")'%(mode,frac,mode,frac))
        exec('l_%s_%s = %s_%s.readlines()'%(mode,frac,mode,frac))

for mode in modes:
    for frac in fracs:
        exec('l0 = l_%s_%s'%(mode,frac))
        exec('l1 = lFix_%s_%s'%(mode,frac))
        for L in l0:
            L = L.replace('\n','')
            l = L.split(' ')
            l = map(float,l)
            if  frac=='nhadron' or frac=='chadron' or frac=='pi0' or frac == 'e' or frac == 'mu' or frac == 'tau' or frac == 'hadron':
                if l[3]: l1.write('%.8g %.8g %.8g %.8g %.8g %.8g'%(l[0],l[1],l[2],l[3],l[4]/l[3],l[5]))
                else: l1.write(L)
                l1.write('\n')
            if frac == 'all':
                if l[3]: l1.write('%.8g %.8g %.8g %.8g %.8g %.8g %.8g'%(l[0],l[1],l[2],l[4],l[3],l[5]/l[3],l[6]))
                else: l1.write(L)
                l1.write('\n')
            if frac == 'other':
                if l[3]: l1.write('%.8g %.8g %.8g %.8g %.8g %.8g'%(l[0],l[1],l[2],l[3],l[4]/l[3],l[5]))
                else: l1.write(L)
                l1.write('\n')

for mode in modes:
    for frac in fracs: exec('lFix_%s_%s.close()'%(mode,frac))
