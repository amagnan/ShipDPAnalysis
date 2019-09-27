from array import array
import os,sys,getopt
import math
lepto=0 
modes = ['meson_pi0','meson_eta','meson_omega','meson_eta1']
fracs = ['e','mu','tau','pi0','nhadron','chadron','hadron','all','other','sum']
#tots = ['other','sum']

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:")

except getopt.GetoptError:
    print 'decay to all SM particles'

for o,a in opts:
    if o in ('-d',): date = a

pathW = "../data/"+date+"/"

f0 = open("mr.txt","r")#dosya lazim
l0 = f0.readlines()

for frac in fracs: exec('l_%s  = open(pathW+"meson_%s.dat","w")'%(frac,frac))

for mode in modes:
    for frac in fracs:
        exec('%s_%s = open(pathW+"%s_%s.dat","r")'%(mode,frac,mode,frac))
        exec('l_%s_%s = %s_%s.readlines()'%(mode,frac,mode,frac))

def find(lines,mass,eps):
    for i in lines:
        k = i.replace('\n','')
        k = k.split(' ')
        if abs(math.log10(float(k[1])) - math.log10(eps)) <0.01 and abs(float(k[0]) - mass)<0.0001: return True
    return False

def find_ratios(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.0001: return float(i[2]), float(i[3]), float(i[4])
    return 0 

def find_allratios(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.0001: return float(i[3]), float(i[4]), float(i[5]), float(i[2])
    return 0

def find_dau(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.0001: return float(i[2]), float(i[3]), float(i[4])
    return 0

def find_N(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.0001: return float(i[2]), float(i[3]), float(i[4]), float(i[5]), float(i[6])
    return 0

def looping(mode,frac,l0,l1):
    l0 = l0.replace('\n','')
    x = l0.split(' ')
    if frac=='e' or frac=='mu' or frac=='tau' or frac=='hadron':
        R  = find_ratios(l1, float(x[0]), float(x[1]))
        exec('N  = find_N(l_%s_sum, float(x[0]), float(x[1]))'%(mode))
        if R and N:
            #print N
            BRn   = R[0]*N[2]
            VPbr  = R[1]*BRn
            GAvp  = R[2]*VPbr
            return BRn, VPbr, GAvp, N[2]
    if frac=='all':
        R  = find_allratios(l1, float(x[0]), float(x[1]))
        #exec('D  = find_dau(%s_other, float(x[0]), float(x[1]))'%(mode)) 
        exec('N  = find_N(l_%s_sum, float(x[0]), float(x[1]))'%(mode))
        if R and N:
            #DAUo  = R[3]*D[0]
            BRn   = R[0]*N[2]
            VPbr  = R[1]*BRn
            GAvp  = R[2]*VPbr
            return BRn, VPbr, GAvp, N[2], N[1]
    if frac=='sum':
        R = find_N(l1, float(x[0]), float(x[1]))
        if R: return R
    if frac=='other':
        R = find_dau(l1, float(x[0]), float(x[1]))
        if R: return R
    return 0

#other ve sum kaldi
#all'un dau'su kaldi.
BR,VP,GA,Dau = 0., 0., 0., 0.
for l in l0:
    l00=l
    l = l.replace('\n','')
    l = l.split(' ')
    mass, eps =float(l[0]), float(l[1])
    for frac in fracs:
        R, V, G, Nr, D = 0., 0., 0., 0., 0.
        NR, DPn, Pnr, Vn, Gn = 0., 0., 0., 0., 0.
        Lp, Lv, Lg = 0., 0., 0.
        for mode in modes:
            fl = 0
            #print mode
            exec('r = looping(mode,frac,l00,l_%s_%s)'%(mode,frac))
            if r:
                exec('k=l_%s'%(frac))
                if  frac=='nhadron' or frac=='chadron' or frac=='pi0' or frac == 'e' or frac == 'mu' or frac == 'tau' or frac == 'hadron':
                    fl = 1
                    R += r[0]
                    V += r[1]
                    G += r[2]
                    Nr += r[3]
                if frac == 'all':
                    fl = 2
                    R += r[0]
                    V += r[1]
                    G += r[2]
                    Nr += r[3]
                    D += r[4]
                if frac == 'sum':
                    fl = -1
                    NR  += r[0]
                    DPn += r[1]
                    Pnr += r[2]
                    Vn  += r[3]
                    Gn  += r[4]
                if frac == 'other':
                    fl = -2
                    Lp += r[0]
                    Lv += r[1]
                    Lg += r[2]
        #print mass, eps, mode, frac
        if fl == 2:
            BR = R/Nr
            if R: VP = V/R
            if V: GA = G/V
            if D: Dau = Nr/D
            k.write("%.8g %.8g %.8g %.8g %.8g %.8g"%(mass,eps,Dau,BR,VP,GA))
            k.write("\n")
        if fl == 1:
            BR = R/Nr
            if R: VP = V/R
            #print frac,mode
            if V: GA = G/V
            k.write("%.8g %.8g %.8g %.8g %.8g"%(mass,eps,BR,VP,GA))
            k.write("\n")
        if fl == -1:
            k.write("%.8g %.8g %.8g %.8g %.8g %.8g %.8g"%(mass,eps,NR,DPn,Pnr,Vn,Gn))
            k.write("\n")
        if fl == -2:
            k.write("%.8g %.8g %.8g %.8g %.8g"%(mass,eps,Lp,Lv,Lg))
            k.write("\n")

for frac in fracs: exec('l_%s.close()'%(frac))
