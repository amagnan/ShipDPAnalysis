from array import array
import os,sys,getopt
import math
from decimal import Decimal

lepto=0 
modes = ['meson','pbrem','qcd']
Fracs = ['e','mu','tau','charged','neutral','all','sum']
fracs = ['e','mu','tau','charged','neutral','all']
ratios = ['Ves','Rec','Geo']

try:
                opts, args = getopt.getopt(sys.argv[1:], "d:")

except getopt.GetoptError:
                print 'decay to all SM particles'

for o,a in opts:
                if o in ('-d',): date = a

pathW = "../data/"+date+"/"

f0 = open("BR_mass.txt","r")#dosya lazim
l0 = f0.readlines()

l_BR    = open(pathW+"BR.dat","w")
l_BRvis = open(pathW+"BR_vis.dat","w")

for mode in modes:
        for frac in Fracs:
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
                                                return BRn, PGn
                if frac=='all':
                                R  = find_allratios(l1, float(x[0]), float(x[1]))
                                #exec('D  = find_dau(%s_other, float(x[0]), float(x[1]))'%(mode)) 
                                exec('N  = find_N(l_%s_sum, float(x[0]), float(x[1]))'%(mode))
                                if R and N:
                                                BRn   = R[1]*N[2]
                                                PGn   = R[0]*N[2]
                                                return BRn, PGn, N[3]
                if frac=='sum':
                                R = find_N(l1, float(x[0]), float(x[1]))
                                if R: return R
                return 0

BR,Dau,PUR,PURN = 0., 0., 0., 0.
for l in l0:
        l00=l
        l = l.replace('\n','')
        l = l.split(' ')
        mass, eps =float(l[0]), float(l[1])
        l_BR.write("%.8g"%mass)
        l_BRvis.write("%.8g"%mass)
        for frac in fracs:
                fl = 0
                R, Pg, PgNr, Nr = 0., 0., 0., 0.
                for mode in modes:
                        exec('N  = find_N(l_%s_sum, mass, eps)'%(mode))
                        if N: Nr += N[2]
                        exec('r = looping(mode,frac,l00,l_%s_%s)'%(mode,frac))
                        if r!=0:
                                #exec('k=l_%s'%(frac))
                                if      frac=='neutral' or frac=='charged' or frac == 'e' or frac == 'mu' or frac == 'tau':
                                        fl = 1
                                        R += r[0]
                                        Pg += r[1]
                                if frac == 'all':
                                        fl = 2
                                        R += r[0]
                                        Pg += r[1]
                                        PgNr += r[2]
                if fl == 2:
                        if Nr!=0.:
                                BR = R/Nr
                                PUR = Pg/Nr
                                PURN = PgNr/Nr
                if fl == 1:
                        if Nr!=0.:
                                BR = R/Nr
                                PUR=Pg/Nr
                if fl !=0:
                        l_BR.write(" %.8g"%BR)
                        l_BRvis.write(" %.8g"%PUR)
                if fl ==0:
                        l_BR.write(" 0.0")
                        l_BRvis.write(" 0.0")
        l_BR.write("\n")
        l_BRvis.write("\n")
l_BR.close()
l_BRvis.close()

"""e_BR,    mu_BR,  tau_BR,         charged_BR,             neutral_BR,             all_BR          = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
e_Vis,  mu_Vis, tau_Vis,        charged_Vis,    neutral_Vis,    all_Vis         = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )

e_eps,  mu_eps, tau_eps,        charged_eps,    neutral_eps,    all_eps         = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
e_mass, mu_mass,tau_mass,       charged_mass,   neutral_mass,   all_mass        = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )

import ROOT

nt=ROOT.TFile(pathW+"Plots.root","RECREATE")
c = ROOT.TCanvas('c', '',1920,1080)

aa=0
for mode in modes:
        aa+=1
        exec('l1=l_%s_all'%mode)
        m=l1[len(l1)-1].replace('\n','')
        #m=l1[len(l1)-1].split(' ')
        #print mode, Decimal(m[0]), Decimal(m[1]), Decimal(m[5]), Decimal(m[6])
        exec('Ves_%s  = ROOT.TH2F("Ves_"+mode, "", 100, 0.,Decimal(float(m[0])+0.2),100,10**-9,10**-4 )'%mode)
        exec('Rec_%s  = ROOT.TH3F("Rec_"+mode, "", 100, 0.,Decimal(float(m[0])+0.2),100,10**-9,10**-4 )'%mode)
        exec('Geo_%s  = ROOT.TH2F("Geo_"+mode, "", 100, 0.,Decimal(float(m[0])+0.2),100,10**-9,10**-4 )'%mode)
        for l in l1:
                l = l.replace('\n','')
                l = l.split(' ')
                geo=Decimal(l[5])*Decimal(l[6])
                #print mode,l[0], Decimal(l[6])
                exec('Ves_%s.Fill(Decimal(l[0]),Decimal(l[1]),Decimal(l[5]))'%mode)
                exec('Rec_%s.Fill(Decimal(l[0]),Decimal(l[1]),Decimal(l[6]))'%mode)
                exec('Geo_%s.Fill(Decimal(l[0]),Decimal(l[1]),geo)'%mode)
        exec('Ves_%s.SetTitle(";m_{#gamma^{D}}(MeV/c^{2});#varepsilon;Vessel Probability")'%mode)
        exec('Ves_%s.Draw("colz")'%mode)
        c.SetLogy()
        c.SetLogz()
        ROOT.gStyle.SetOptStat(0000)
        c.Modified()
        c.Update()
        c.Write()
        c.Print(pathW+"Ves_"+mode+".pdf")
        exec('Rec_%s.SetTitle(";m_{#gamma^{D}}(MeV/c^{2});#varepsilon;Reconstruction Efficiency")'%mode)
        exec('Rec_%s.Draw("colz")'%mode)
        c.SetLogy()
        c.SetLogz()
        ROOT.gStyle.SetOptStat(0000)
        c.Modified()
        c.Update()
        c.Write()
        c.Print(pathW+"Rec_"+mode+".pdf")
        exec('Geo_%s.SetTitle(";m_{#gamma^{D}}(MeV/c^{2});#varepsilon;Final Acceptance")'%mode)
        exec('Geo_%s.Draw("colz")'%mode)
        c.SetLogy()
        c.SetLogz()
        ROOT.gStyle.SetOptStat(0000)
        c.Modified()
        c.Update()
        c.Write()
        c.Print(pathW+"Geo_"+mode+".pdf")
nt.Write()
nt.Print()
nt.Close()
print aa
#l_BR   = open(pathW+"BR.dat","r")
#l_BRvis = open(pathW+"BR_vis.dat","r")
#l_BR   = l_BR.readlines()
#l_BRvis = l_BRvis.readlines()"""
