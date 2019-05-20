import ROOT as r
from array import array
import math
import numpy
import os,sys,getopt
try:
    opts, args = getopt.getopt(sys.argv[1:], "d:")

except getopt.GetoptError:
    print 'no file'
    sys.exit()

for o,a in opts:
    if o in ('-d',): date = a

def find_ratios(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.1 and abs(float(i[0]) - mass)<0.0001: return float(i[2]), float(i[3]), float(i[4])
    else: return 0 

def find_allratios(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.1 and abs(float(i[0]) - mass)<0.0001: return float(i[3]), float(i[4]), float(i[5])
    else: return 0

def find_mass(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.1 and abs(float(i[0]) - mass)<0.0001: return float(i[0])
    else: return 0


def find_dau(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.1 and abs(float(i[0]) - mass)<0.0001: return float(i[2])
    else: return 0

def find_N(lines,mass,eps):
    for i in lines:
        i = i.replace('\n','')
        i = i.split(' ')
        if abs(math.log10(float(i[1])) - math.log10(eps)) <0.1 and abs(float(i[0]) - mass)<0.0001: return float(i[4])
    else: return 0

def looping(line,itsmass,itseps,itsbr,itsves,itsreco):
    for i in line:
        i = i.replace('\n','')
        i = i.split(' ')
        #print i[1],i[4]
        itsmass.append(float(i[0]))
        itseps.append(float(i[1]))
        itsbr.append(float(i[2]))
        itsves.append(float(i[3]))
        itsreco.append(float(i[4]))

def looping_all(line,itsmass,itseps,itsdau,itsbr,itsvtx,itsreco):
    for i in line:
        i = i.replace('\n','')
        i = i.split(' ')      
        #print i[0],i[1]
        itsmass.append(float(i[0]))
        itseps.append(float(i[1]))
        itsdau.append(float(i[2]))
        itsbr.append(float(i[3]))
        itsreco.append(float(i[5]))
        itsvtx.append(float(i[4]))
 
def looping_oth(line,itsmass,itseps,itsdau):
     for i in line:
        i = i.replace('\n','')
        i = i.split(' ') 
        itsmass.append(float(i[0]))
        itseps.append(float(i[1]))
        itsdau.append(float(i[2]))

def looping_N(line,itsmass,itseps,itssum):
     for i in line:
        i = i.replace('\n','')
        i = i.split(' ') 
        itsmass.append(float(i[0]))
        itseps.append(float(i[1]))
        itssum.append(float(i[4]))

def combined_mode(mode,eps,l_meson_mode,l_pbrem_mode,l_qcd_mode,combined_mode_mass,combined_mode_br):

    for x in l_meson_mode:
        x = x.replace('\n','')
        x = x.split(' ')
        if mode!='all': R1 = find_ratios(l_meson_mode,float(x[0]),eps)
        if mode=='all': R1 = find_allratios(l_meson_mode,float(x[0]),eps)
        y = find_mass(l_pbrem_mode,float(float(x[0])),float(eps))
        if y:
            if y!=float(x[0]): print 'something is wrong!'
            N1 = find_N(l_meson_sum,float(x[0]),eps)
            N2 = find_N(l_pbrem_sum,float(x[0]),eps)
            if mode!='all': R2 = find_ratios(l_pbrem_mode,float(x[0]),eps)
            if mode=='all': R2 = find_allratios(l_pbrem_mode,float(x[0]),eps)
            BR = (N1*R1[0]+N2*R2[0])/(N1+N2)
        if not y:
            BR = R1[0]
        combined_mode_mass.append(float(x[0]))
        combined_mode_br.append(BR)

    for k in l_pbrem_mode:
        k = k.replace('\n','')
        k = k.split(' ')
        if not float(k[0]) in combined_mode_mass:
            if mode!='all': R2 = find_ratios(l_pbrem_mode,float(k[0]),eps)
            if mode=='all': R2 = find_allratios(l_pbrem_mode,float(k[0]),eps)
            z = find_mass(l_qcd_mode,float(float(k[0])),float(eps))
            if z:
                if z!=float(k[0]): print 'something is wrong!'
                N2 = find_N(l_pbrem_sum,float(k[0]),eps)
                N3 = find_N(l_qcd_sum,float(k[0]),eps)
                if mode!='all': R3 = find_ratios(l_qcd_mode,float(k[0]),eps)
                if mode=='all': R3 = find_allratios(l_qcd_mode,float(k[0]),eps)
                BR = (N2*R2[0]+N3*R3[0])/(N2+N3)
            if not z:
                BR = R2[0]
            combined_mode_mass.append(float(k[0]))
            combined_mode_br.append(BR)

    for l in l_qcd_mode:
        l = l.replace('\n','')
        l = l.split(' ')
        if not float(l[0]) in combined_mode_mass:
            if mode!='all': R3 = find_ratios(l_qcd_mode,float(l[0]),eps)
            if mode=='all': R3 = find_allratios(l_qcd_mode,float(l[0]),eps)
            BR = R3[0]
            combined_mode_mass.append(float(l[0]))
            combined_mode_br.append(BR)
def Sort(mass,eps):
    m=mass
    e=eps
    mass, eps = array('d'), array('d')
    s=numpy.argsort(m)
    for a in s:
        mass.append(float(m[a]))
        eps.append(float(e[a]))

def plots(prod,e_mass,mu_mass,tau_mass,hadron_mass,all_mass,e_br,mu_br,tau_br,hadron_br,all_br):

    br  = r.TMultiGraph()
    
    if len(e_mass)!=0.:
        Sort(e_mass,e_br)
        #print e_mass
        n=len(e_mass)
        #print e_br
        e_Br=r.TGraph(n,e_mass,e_br)
        e_Br.SetName("e_Br")
        e_Br.SetMarkerColor(2)
        e_Br.SetMarkerStyle(20)
        e_Br.SetLineColor(2)
        e_Br.SetLineWidth(2)
        e_Br.SetTitle("e^{+} + e^{-}") 
        e_Br.SetFillStyle(0) 
        e_Br.SetDrawOption("AP")
        br.Add(e_Br)
    
    if len(mu_mass)!=0.:
        Sort(mu_mass,mu_br)
        mu_Br=r.TGraph(len(mu_mass),mu_mass,mu_br)
        mu_Br.SetName("mu_Br")
        mu_Br.SetMarkerColor(4)
        mu_Br.SetMarkerStyle(21)
        mu_Br.SetLineColor(4)
        mu_Br.SetLineWidth(2)
        mu_Br.SetTitle("#mu^{+} + #mu^{-}") 
        mu_Br.SetFillStyle(0) 
        mu_Br.SetDrawOption("AP")
        br.Add(mu_Br)
    
    if len(tau_mass)!=0.:
        Sort(tau_mass,tau_br)
        tau_Br=r.TGraph(len(tau_mass),tau_mass,tau_br)
        tau_Br.SetName("tau_Br")
        tau_Br.SetMarkerColor(3)
        tau_Br.SetMarkerStyle(22)
        tau_Br.SetLineColor(3)
        tau_Br.SetLineWidth(2)
        tau_Br.SetTitle("#tau^{+} + #tau^{-}") 
        tau_Br.SetFillStyle(0) 
        tau_Br.SetDrawOption("AP")
        br.Add(tau_Br)
    
    if len(hadron_mass)!=0.:
        Sort(hadron_mass,hadron_br)
        hadron_Br=r.TGraph(len(hadron_mass),hadron_mass,hadron_br)
        #print hadron_mass
        hadron_Br.SetName("hadron_Br")
        hadron_Br.SetMarkerColor(6)
        hadron_Br.SetMarkerStyle(34)
        hadron_Br.SetLineColor(6)
        hadron_Br.SetLineWidth(2)
        hadron_Br.SetTitle("hadrons") 
        hadron_Br.SetFillStyle(0) 
        hadron_Br.SetDrawOption("AP")
        br.Add(hadron_Br)
    
    if len(all_mass)!=0.:
        Sort(all_mass,all_br)
        all_Br=r.TGraph(len(all_mass),all_mass,all_br)
        all_Br.SetName("all_Br")
        all_Br.SetMarkerColor(1)
        all_Br.SetMarkerStyle(29)
        all_Br.SetLineColor(1)
        all_Br.SetLineWidth(2)
        all_Br.SetTitle("All modes") 
        all_Br.SetFillStyle(0) 
        all_Br.SetDrawOption("AP")
        br.Add(all_Br)

    br.SetTitle(";Mass(GeV);Branching Ratios")
    br.Draw("AP")
    #xaxis=br.GetXaxis()
    #xaxis.SetLimits(0.,6.)
    c1.BuildLegend()
    c1.SetLogx()
    c1.Modified()
    c1.Update()
    c1.Write()
    c1.Print(pathW+prod+"_br.pdf")

pathR = "../data/"+date+"/"
prods=['meson','pbrem','qcd']
pathW="../plots/"+date+"/"
modes=['e','mu','tau','hadron','all']
 
hfile   =   r.TFile(pathW+"plots.root","RECREATE","")
c1      =   r.TCanvas('c1', '',1920,1080)


parameters=['mass','eps']
variables=['br','vessel','reco']
#other,sum/N ve dau elle verildi

for n,mode in enumerate(modes):
    exec('combined_%s_mass = array( "d" )' %mode)
    exec('combined_%s_br = array( "d" )' %mode)

for i,prod in enumerate(prods):
    #print i,prod#istisnalar bitti
    exec('%s_other_dau = array( "d" )' %prod)
    exec('%s_all_dau = array( "d" )' %prod)
    exec('%s_sum_N = array( "d" )' %prod)

    exec('f_%s_other = open(pathR+prod+"_Ana_other.dat","r")' %prod)
    exec('f_%s_sum = open(pathR+prod+"_Ana_sum.dat","r")' %prod)

    exec('l_%s_other = f_%s_other.readlines()' %(prod,prod))
    exec('l_%s_sum = f_%s_sum.readlines()' %(prod,prod))

    exec('t_%s = r.TTree("%s","")'%(prod,prod))

    for m,parameter in enumerate(parameters):
        #print m,parameter
        exec('%s_other_%s = array( "d" )' %(prod,parameter))
        exec('%s_sum_%s = array( "d" )' %(prod,parameter))


    for n,mode in enumerate(modes):#dosyalar acildi,line by line okundu
        #print n,mode
        exec('f_%s_%s = open(pathR+prod+"_Ana_"+mode+".dat","r")' %(prod,mode))
        exec('l_%s_%s = f_%s_%s.readlines()' %(prod,mode,prod,mode))


        for m,parameter in enumerate(parameters):
            #print m,parameter
            exec('%s_%s_%s = array( "d" )' %(prod,mode,parameter))

        for k,variable in enumerate(variables):
            #print k,variable
            exec('%s_%s_%s = array( "d" )' %(prod,mode,variable))
    
    exec('looping_all(l_%s_all,%s_all_mass,%s_all_eps,%s_all_dau,%s_all_br,%s_all_vessel,%s_all_reco)' %(prod,prod,prod,prod,prod,prod,prod))
    exec('looping_oth(l_%s_other,%s_other_mass,%s_other_eps,%s_other_dau)' %(prod,prod,prod,prod))
    exec('looping_N(l_%s_sum,%s_sum_mass,%s_sum_eps,%s_sum_N)' %(prod,prod,prod,prod))
    exec('plots(prod,%s_e_mass,%s_mu_mass,%s_tau_mass,%s_hadron_mass,%s_all_mass,%s_e_br,%s_mu_br,%s_tau_br,%s_hadron_br,%s_all_br)'%(prod,prod,prod,prod,prod,prod,prod,prod,prod,prod))
    for n,mode in enumerate(modes):
        if n!=4:
            #print n,mode,'for looping'
            exec('looping(l_%s_%s,%s_%s_mass,%s_%s_eps,%s_%s_br,%s_%s_vessel,%s_%s_reco)' %(prod,mode,prod,mode,prod,mode,prod,mode,prod,mode,prod,mode))
        else: continue


    for n,mode in enumerate(modes):
        #print n,mode,"h2ler icin"
        for k,variable in enumerate(variables):
            #print k,variable,"h2ler icin"
            #exec('if len(%s_%s_%s) > 0:'%(prod,mode,variable))
            #    if variable == "br":
            exec('if len(%s_%s_%s) > 0 and variable == "br": h2_%s_%s_%s = r.TH2F( "h2_%s_%s_%s", ";log(mass(GeV));log(#varepsilon);Branching Ratios of %s", int(math.sqrt(len(%s_%s_mass))),  math.log10(min(%s_%s_mass)), math.log10(max(%s_%s_mass)), int(math.sqrt(len(%s_%s_eps))), math.log10(min(%s_%s_eps)), math.log10(max(%s_%s_eps)))'%(prod,mode,variable,prod,mode,variable,prod,mode,variable,mode,prod,mode,prod,mode,prod,mode,prod,mode,prod,mode,prod,mode))
            #   if variable == "vessel":
            exec('if len(%s_%s_%s) > 0 and variable == "vessel": h2_%s_%s_%s = r.TH2F( "h2_%s_%s_%s", ";log(mass(GeV));log(#varepsilon);Vessel Probability of %s", int(math.sqrt(len(%s_%s_mass))),  math.log10(min(%s_%s_mass)), math.log10(max(%s_%s_mass)), int(math.sqrt(len(%s_%s_eps))), math.log10(min(%s_%s_eps)), math.log10(max(%s_%s_eps)))'%(prod,mode,variable,prod,mode,variable,prod,mode,variable,mode,prod,mode,prod,mode,prod,mode,prod,mode,prod,mode,prod,mode))
            #    if variable == "reco":
            exec('if len(%s_%s_%s) > 0 and variable == "reco": h2_%s_%s_%s = r.TH2F( "h2_%s_%s_%s", ";log(mass(GeV));log(#varepsilon);Reconstruction Efficiency of %s", int(math.sqrt(len(%s_%s_mass))),  math.log10(min(%s_%s_mass)), math.log10(max(%s_%s_mass)), int(math.sqrt(len(%s_%s_eps))), math.log10(min(%s_%s_eps)), math.log10(max(%s_%s_eps)))'%(prod,mode,variable,prod,mode,variable,prod,mode,variable,mode,prod,mode,prod,mode,prod,mode,prod,mode,prod,mode,prod,mode))

    exec('if len(%s_other_dau) > 0: h2_%s_other_dau = r.TH2F( "h2_%s_other_dau", ";log(mass(GeV));log(#varepsilon);Not proper Daughter Events/All DP Events", int(math.sqrt(len(%s_other_mass))),  math.log10(min(%s_other_mass)), math.log10(max(%s_other_mass)), int(math.sqrt(len(%s_other_eps))), math.log10(min(%s_other_eps)), math.log10(max(%s_other_eps)))'%(prod,prod,prod,prod,prod,prod,prod,prod,prod))
    exec('if len(%s_all_dau) > 0: h2_%s_all_dau = r.TH2F( "h2_%s_all_dau", ";log(mass(GeV));log(#varepsilon);Proper Daughter Events/All DP Events", int(math.sqrt(len(%s_all_mass))),  math.log10(min(%s_all_mass)), math.log10(max(%s_all_mass)), int(math.sqrt(len(%s_all_eps))), math.log10(min(%s_all_eps)), math.log10(max(%s_all_eps)))'%(prod,prod,prod,prod,prod,prod,prod,prod,prod))

    exec('nFile=len(%s_all_mass)'%prod)
    #vectorler ve branchler burada tanimlanacak
    exec('v_%s_all_dau = r.std.vector("double")()'%prod)
    exec('v_%s_other_dau = r.std.vector("double")()'%prod)
    exec('v_%s_sum_N = r.std.vector("double")()'%prod)
    exec('t_%s.Branch("dau_all_"+prod, v_%s_all_dau)' %(prod,prod))
    exec('t_%s.Branch("dau_other_"+prod, v_%s_other_dau)'%(prod,prod))
    exec('t_%s.Branch("N_sum_"+prod, v_%s_sum_N)'%(prod,prod))
    exec('v2_%s_all_dau = r.std.vector(r.TVector3)()'%prod)
    exec('v2_%s_other_dau = r.std.vector(r.TVector3)()'%prod)
    exec('t_%s.Branch("Twodau_all_"+prod, v2_%s_all_dau)'%(prod,prod))
    exec('t_%s.Branch("Twodau_other_"+prod, v2_%s_other_dau)'%(prod,prod))
    for m,parameter in enumerate(parameters):
        exec('v_%s_sum_%s = r.std.vector("double")()'%(prod,parameter))
        exec('v_%s_other_%s = r.std.vector("double")()'%(prod,parameter))
        exec('t_%s.Branch(parameter+"_other_"+prod, v_%s_other_%s)'%(prod,prod,parameter))
        exec('t_%s.Branch(parameter+"_sum_"+prod, v_%s_sum_%s)'%(prod,prod,parameter))
    for n,mode in enumerate(modes):
        for k,variable in enumerate(variables):
            exec('v_%s_%s_%s = r.std.vector("double")()'%(prod,mode,variable))
            exec('t_%s.Branch(variable+"_"+mode+"_"+prod, v_%s_%s_%s)'%(prod,prod,mode,variable))
            exec('v2_%s_%s_%s = r.std.vector(r.TVector3)()'%(prod,mode,variable))
            exec('t_%s.Branch("Two"+variable+"_"+mode+"_"+prod, v2_%s_%s_%s)'%(prod,prod,mode,variable))
        for m,parameter in enumerate(parameters):
            exec('v_%s_%s_%s = r.std.vector("double")()'%(prod,mode,parameter))
            exec('t_%s.Branch(parameter+"_"+mode+"_"+prod, v_%s_%s_%s)'%(prod,prod,mode,parameter))

    for j in range(nFile):
        exec('v_%s_all_dau.clear()'%prod)
        exec('v_%s_other_dau.clear()'%prod)
        exec('v2_%s_all_dau.clear()'%prod)
        exec('v2_%s_other_dau.clear()'%prod)
        exec('v_%s_sum_N.clear()'%prod)
        for m,parameter in enumerate(parameters):
            exec('v_%s_sum_%s.clear()'%(prod,parameter))
            exec('v_%s_other_%s.clear()'%(prod,parameter))
        for n,mode in enumerate(modes):
            for k,variable in enumerate(variables):
                exec('v_%s_%s_%s.clear()'%(prod,mode,variable))
                exec('v2_%s_%s_%s.clear()'%(prod,mode,variable))
            for m,parameter in enumerate(parameters):
                exec('v_%s_%s_%s.clear()'%(prod,mode,parameter))
        #exec('%s_all = %s_all_br[j],%s_all_vessel[j],%s_all_reco[j]'%(prod,prod,prod,prod))
        #exec('print %s_all[0],%s_all[1],%s_all[2]'%(prod,prod,prod))
        #exec('if %s_all_br[j]>1.0: print %s_all_mass[j], %s_all_eps[j], %s_all_br[j]' %(prod,prod,prod,prod)) 
        exec('mass, eps = %s_all_mass[j], %s_all_eps[j]'%(prod,prod)) 
        for n,mode in enumerate(modes):
            exec('%s_%s = find_ratios(l_%s_%s, mass, eps)' %(prod,mode,prod,mode) )
        exec('%s_other = find_dau(l_%s_other, mass, eps)' %(prod,prod) )
        exec('%s_sum = find_N(l_%s_sum, mass, eps)' %(prod,prod) )

        exec('if %s_other: v_%s_other_dau.push_back(%s_other)'%(prod,prod,prod))
        exec('if %s_sum: v_%s_sum_N.push_back(%s_sum)'%(prod,prod,prod))
        exec('if %s_other: v2_%s_other_dau.push_back(r.TVector3(math.log10(mass),math.log10(eps),%s_other))'%(prod,prod,prod))
        exec('if %s_other: h2_%s_other_dau.Fill(math.log10(mass),math.log10(eps),%s_other)'%(prod,prod,prod))
        for m,parameter in enumerate(parameters):
            exec('if %s_sum: v_%s_sum_%s.push_back(math.log10(%s))'%(prod,prod,parameter,parameter))
            exec('if %s_other: v_%s_other_%s.push_back(math.log10(%s))'%(prod,prod,parameter,parameter))
        for n,mode in enumerate(modes):
            if mode =="all": continue
            for k,variable in enumerate(variables):
                #print mode,k,variable
                exec('if %s_%s: v_%s_%s_%s.push_back(%s_%s[k])'%(prod,mode,prod,mode,variable,prod,mode))
                exec('if %s_%s: v2_%s_%s_%s.push_back(r.TVector3(math.log10(mass),math.log10(eps),%s_%s[k]))'%(prod,mode,prod,mode,variable,prod,mode))
                exec('if %s_%s: h2_%s_%s_%s.Fill(math.log10(mass),math.log10(eps),%s_%s[k])'%(prod,mode,prod,mode,variable,prod,mode))
            for m,parameter in enumerate(parameters):
                exec('if %s_%s: v_%s_%s_%s.push_back(math.log10(%s))'%(prod,mode,prod,mode,parameter,parameter))
        exec('v_%s_all_dau.push_back(%s_all_dau[j])'%(prod,prod))
        exec('v2_%s_all_dau.push_back(r.TVector3(math.log10(mass),math.log10(eps),%s_all_dau[j]))'%(prod,prod))
        exec('if %s_all_dau[j]: h2_%s_all_dau.Fill(math.log10(mass),math.log10(eps),%s_all_dau[j])'%(prod,prod,prod))
        for k,variable in enumerate(variables):
            exec('if %s_all_%s: v_%s_all_%s.push_back(%s_all_%s[j])'%(prod,variable,prod,variable,prod,variable))
            exec('if %s_all_%s: v2_%s_all_%s.push_back(r.TVector3(math.log10(mass),math.log10(eps),%s_all_%s[j]))'%(prod,variable,prod,variable,prod,variable))
            exec('if %s_all_%s: h2_%s_all_%s.Fill(math.log10(mass),math.log10(eps),%s_all_%s[j])'%(prod,variable,prod,variable,prod,variable))
        for m,parameter in enumerate(parameters):
             exec('v_%s_all_%s.push_back(math.log10(%s))'%(prod,parameter,parameter))
        #exec('if %s_all_br[j]!=1.0: print %s_all_mass[j], %s_all_eps[j],%s_all_br[j]' %(prod,prod,prod,prod))
        exec('t_%s.Fill()'%prod)

    exec('if %s_other: h2_%s_other_dau.SetAxisRange(min(%s_other_dau),max(%s_other_dau),"Z")' %(prod,prod,prod,prod))
    exec('if %s_other: h2_%s_other_dau.Draw("colz")' %(prod,prod))
    #exec('if %s_other: h2_%s_other_dau.Write("colz")' %(prod,prod))
    r.gStyle.SetOptStat(0000)
    c1.Modified()
    c1.Update()
    c1.Write()
    exec('c1.Print(pathW+"%s_other_dau.pdf")'%prod)
    for n,mode in enumerate(modes):
        for k,variable in enumerate(variables):
            exec('if %s_%s: h2_%s_%s_%s.SetAxisRange(min(%s_%s_%s),max(%s_%s_%s),"Z")' %(prod,mode,prod,mode,variable,prod,mode,variable,prod,mode,variable))
            exec('if %s_%s: h2_%s_%s_%s.Draw("colz")' %(prod,mode,prod,mode,variable))
            #exec('if %s_%s: h2_%s_%s_%s.Write()' %(prod,mode,prod,mode,variable))
            r.gStyle.SetOptStat(0000)
            c1.Modified()
            c1.Update()
            c1.Write()
            exec('c1.Print(pathW+"%s_%s_%s.pdf")'%(prod,mode,variable))
    exec('if %s_all_dau: h2_%s_all_dau.SetAxisRange(min(%s_all_dau),max(%s_all_dau),"Z")' %(prod,prod,prod,prod))
    exec('if %s_all_dau: h2_%s_all_dau.Draw("colz")' %(prod,prod))
    #exec('if %s_all_dau: h2_%s_all_dau.Write()' %(prod,prod))
    r.gStyle.SetOptStat(0000)
    c1.Modified()
    c1.Update()
    c1.Write()
    exec('c1.Print(pathW+"%s_all_dau.pdf")'%prod)
for n,mode in enumerate(modes):
    exec('combined_mode(mode,1e-06,l_meson_%s,l_pbrem_%s,l_qcd_%s,combined_%s_mass,combined_%s_br)'%(mode,mode,mode,mode,mode))
plots('combined',combined_e_mass,combined_mu_mass,combined_tau_mass,combined_hadron_mass,combined_all_mass,combined_e_br,combined_mu_br,combined_tau_br,combined_hadron_br,combined_all_br)
hfile.Write()
hfile.Close()
