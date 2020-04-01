from array import array
import os,sys,getopt
import math
lepto=0 
modes = ['meson_pi0','meson_eta','meson_omega','meson_eta1','meson_eta11']
#modes = ['meson_eta','meson_omega','meson_eta1']
fracs = ['e','mu','tau','neutral','charged','all','other','sum']
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
        if abs(math.log10(float(k[1])) - math.log10(eps)) <0.01 and abs(float(k[0]) - mass)<0.00001: return True
    return False

def find_ratios(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.00001: return float(i[2]), float(i[3]), float(i[4]), float(i[5])
    return 0 

def find_allratios(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.00001: return float(i[3]), float(i[4]), float(i[5]), float(i[6]), float(i[2])
    return 0

def find_dau(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.00001: return float(i[2]), float(i[3]), float(i[4]), float(i[5])
    return 0

def find_N(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.01 and abs(float(i[0]) - mass)<0.00001: return float(i[2]), float(i[3]), float(i[4]), float(i[5]), float(i[6]), float(i[7])
    return 0

def looping(mode,frac,l0,l1):
    l0 = l0.replace('\n','')
    x = l0.split(' ')
    if frac=='e' or frac=='mu' or frac=='tau' or frac=='neutral' or frac=='charged':
        R  = find_ratios(l1, float(x[0]), float(x[1]))
        exec('N  = find_N(l_%s_sum, float(x[0]), float(x[1]))'%(mode))
        if R and N:
            #print N
            BRn   = R[0]*N[2]
            PGn   = R[1]*N[2]
            VPbr  = R[2]*PGn
            GAvp  = R[3]*VPbr
            return BRn, PGn, VPbr, GAvp, N[2]
    if frac=='all':
        R  = find_allratios(l1, float(x[0]), float(x[1]))
        #exec('D  = find_dau(%s_other, float(x[0]), float(x[1]))'%(mode)) 
        exec('N  = find_N(l_%s_sum, float(x[0]), float(x[1]))'%(mode))
        if R and N:
            #DAUo  = R[3]*D[0]
            BRn   = R[1]*N[2]
            PGn   = R[0]*N[2]
            VPbr  = R[2]*PGn
            GAvp  = R[3]*VPbr
            return BRn, PGn, VPbr, GAvp, N[2], N[1], N[3]
    if frac=='sum':
        R = find_N(l1, float(x[0]), float(x[1]))
        if R: return R
    if frac=='other':
        R = find_dau(l1, float(x[0]), float(x[1]))
        if R: return R
    return 0

#other ve sum kaldi
#all'un dau'su kaldi.
BR,VP,GA,Dau,PUR = 0., 0., 0., 0., 0.
for l in l0:
    l00=l
    l = l.replace('\n','')
    l = l.split(' ')
    mass, eps =float(l[0]), float(l[1])
    for frac in fracs:
        R, V, G, Nr, D = 0., 0., 0., 0., 0.
        NR, DPn, Pg, PgNr, Pnr,PurNr, Vn, Gn = 0., 0., 0., 0., 0., 0., 0., 0.
        Lp, Lv, Lg, Lpur = 0., 0., 0., 0.
        for mode in modes:
            fl = 0
            #print mode
            exec('r = looping(mode,frac,l00,l_%s_%s)'%(mode,frac))
            if r:
                exec('k=l_%s'%(frac))
                if  frac=='neutral' or frac=='charged' or frac == 'e' or frac == 'mu' or frac == 'tau':
                    fl = 1
                    R += r[0]
                    V += r[2]
                    G += r[3]
                    Pg += r[1]
                    Nr += r[4]
                if frac == 'all':
                    fl = 2
                    R += r[0]
                    V += r[2]
                    G += r[3]
                    Nr += r[4]
                    D += r[5]
                    Pg += r[1]
                    PgNr += r[6]
                if frac == 'sum':
                    fl = -1
                    NR  += r[0]
                    DPn += r[1]
                    Pnr += r[2]
                    PurNr  += r[3]
                    Vn  += r[4]
                    Gn  += r[5]
                if frac == 'other':
                    fl = -2
                    Lp += r[0]
                    Lpur += r[1]
                    Lv += r[2]
                    Lg += r[3]
        #print mass, eps, mode, frac
        if fl == 2:
            BR = R/Nr
            PUR = Pg/Nr
            if Pg: VP=V/Pg
            #if R:
                #VP = V/R
                #PUR = Pg/R
            if V!=0.: GA = G/V
            if V==0.: GA=0.
            if D: Dau = Nr/D
            if PgNr: PURN = PgNr/Nr
            k.write("%.8g %.8g %.8g %.8g %.8g %.8g %.8g"%(mass,eps,Dau,PURN,BR,VP,GA))
            k.write("\n")
            if abs(PURN - PUR)>0.00001: print PURN, PUR
        if fl == 1:
            BR = R/Nr
            PUR=Pg/Nr
            if Pg: VP=V/Pg
            #if R:
                #PUR = Pg/R
                #VP = V/R
            #print frac,mode
            if V!=0.: GA = G/V
            if V==0.: GA=0
            k.write("%.8g %.8g %.8g %.8g %.8g %.8g"%(mass,eps,BR,PUR,VP,GA))
            k.write("\n")
        if fl == -1:
            k.write("%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g"%(mass,eps,NR,DPn,Pnr,PurNr,Vn,Gn))
            k.write("\n")
        if fl == -2:
            k.write("%.8g %.8g %.8g %.8g %.8g %.8g"%(mass,eps,Lp,Lpur,Lv,Lg))
            k.write("\n")

for frac in fracs: exec('l_%s.close()'%(frac))
