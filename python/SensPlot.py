import ROOT
from array import array
import os,sys,getopt
import math
from decimal import Decimal
lepto=0 
try:
                opts, args = getopt.getopt(sys.argv[1:], "d:")

except getopt.GetoptError:
                print 'decay to all SM particles'

for o,a in opts:
                if o in ('-d',): date = a

pathW = "../data/"+date+"/"
nt=ROOT.TFile(pathW+"Plots.root","RECREATE")
c = ROOT.TCanvas('c', '',1920,1080)

#For the axis titles:
ROOT.gStyle.SetTitleColor(1, "XYZ")
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetTitleSize(0.045, "XYZ")

ROOT.gStyle.SetTitleXOffset(1.2)
ROOT.gStyle.SetTitleYOffset(0.9)

#For the axis labels:
ROOT.gStyle.SetLabelColor(1, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetLabelOffset(0.005, "XYZ")
ROOT.gStyle.SetLabelSize(0.035, "XYZ")

#For the axis:
ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetStripDecimals(True)
ROOT.gStyle.SetTickLength(0.02, "XYZ")
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetPadTickX(1)              #To get tick marks on the opposite side of the frame
ROOT.gStyle.SetPadTickY(1)
#c.SetLeftMargin(0.20)
#c.SetTopMargin(0.15)
c.SetBottomMargin(0.15)
c.SetRightMargin(0.20)
c.SetGrid()
pbrem_all = open(pathW+"pbrem_all.dat","r")
l_pbrem_all = pbrem_all.readlines()

l1=l_pbrem_all
m=l1[len(l1)-1].replace('\n','')
m=l1[len(l1)-1].split(' ')
Ves_pbrem       = ROOT.TH2D("Ves_pbrem", "", 60, -3.,0.7,60,-19.,-6.)
Rec_pbrem       = ROOT.TH2D("Rec_pbrem", "", 60, -3.,0.7,60,-19.,-6.)
Geo_pbrem       = ROOT.TH2D("Geo_pbrem", "", 60, -3.,0.7,60,-19.,-6.)
#k= open(pathW+"pbrem_A.dat","w")
for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        geo=Decimal(float(lm[5])*100.)*Decimal(lm[6])
        if float(lm[5])!=0.: Ves_pbrem.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), Decimal(float(lm[5])*100.))
        if float(lm[6])!=0.: Rec_pbrem.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), Decimal(float(lm[6])*100.))
        if geo!=0.: Geo_pbrem.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), geo)

Ves_pbrem.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Vessel Probability(%)")
Ves_pbrem.SetAxisRange(0.,6.,"Z")
Ves_pbrem.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Ves_pbrem.pdf")

Rec_pbrem.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Reconstruction Efficiency(%)")
Rec_pbrem.SetAxisRange(0.,100.,"Z")
Rec_pbrem.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Rec_pbrem.pdf")

Geo_pbrem.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Final Acceptance(%)")
Geo_pbrem.SetAxisRange(0.,6.,"Z")
Geo_pbrem.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Geo_pbrem.pdf")

meson_all = open(pathW+"meson_all.dat","r")
l_meson_all = meson_all.readlines()

l1=l_meson_all
m=l1[len(l1)-1].replace('\n','')
m=l1[len(l1)-1].split(' ')
Ves_meson       = ROOT.TH2D("Ves_meson", "", 60, -3.,0.1,60,-19.,-6.)
Rec_meson       = ROOT.TH2D("Rec_meson", "", 60, -3.,0.1,60,-19.,-6.)
Geo_meson       = ROOT.TH2D("Geo_meson", "", 60, -3.,0.1,60,-19.,-6.)
#k= open(pathW+"meson_A.dat","w")
for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        geo=Decimal(float(lm[5])*100.)*Decimal(lm[6])
        if float(lm[5])!=0.: Ves_meson.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), Decimal(float(lm[5])*100.))
        if float(lm[6])!=0.: Rec_meson.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), Decimal(float(lm[6])*100.))
        if geo!=0.: Geo_meson.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), geo)

Ves_meson.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Vessel Probability(%)")
Ves_meson.SetAxisRange(0.,6.,"Z")
Ves_meson.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Ves_meson.pdf")

Rec_meson.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Reconstruction Efficiency(%)")
Rec_meson.SetAxisRange(0.,100.,"Z")
Rec_meson.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Rec_meson.pdf")

Geo_meson.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Final Acceptance(%)")
Geo_meson.SetAxisRange(0.,6.,"Z")
Geo_meson.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Geo_meson.pdf")

qcd_all = open(pathW+"qcd_all.dat","r")
l_qcd_all = qcd_all.readlines()

l1=l_qcd_all
m=l1[len(l1)-1].replace('\n','')
m=l1[len(l1)-1].split(' ')
Ves_qcd = ROOT.TH2D("Ves_qcd", "", 60, 0.,1.1,60,-18.,-11.)
Rec_qcd = ROOT.TH2D("Rec_qcd", "", 60, 0.,1.1,60,-18.,-11.)
Geo_qcd = ROOT.TH2D("Geo_qcd", "", 60, 0.,1.1,60,-18.,-11.)
#k= open(pathW+"qcd_A.dat","w")
for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        geo=Decimal(float(lm[5])*100.)*Decimal(lm[6])
        if float(lm[5])!=0.: Ves_qcd.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), Decimal(float(lm[5])*100.))
        if float(lm[6])!=0.: Rec_qcd.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), Decimal(float(lm[6])*100.))
        if geo!=0.: Geo_qcd.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])**2.), geo)

Ves_qcd.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Vessel Probability(%)")
Ves_qcd.SetAxisRange(0.,6.,"Z")
Ves_qcd.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Ves_qcd.pdf")

Rec_qcd.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Reconstruction Efficiency(%)")
Rec_qcd.SetAxisRange(0.,100.,"Z")
Rec_qcd.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Rec_qcd.pdf")

Geo_qcd.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});#varepsilon^{2};Final Acceptance(%)")
Geo_qcd.SetAxisRange(0.,6.,"Z")
Geo_qcd.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Geo_qcd.pdf")

ALPACA_ratio = open("../data/alp/ALPACA_ratio.dat","r")
l_ALPACA_ratio = ALPACA_ratio.readlines()

l1=l_ALPACA_ratio
m=l1[len(l1)-1].replace('\n','')
m=l1[len(l1)-1].split(' ')

Ves_ALPACA = ROOT.TH2D("Ves_ALPACA", "", 45, -2.5,0.4,30,-7.1,-5.5)
Rec_ALPACA = ROOT.TH2D("Rec_ALPACA", "", 45, -2.5,0.4,30,-7.1,-5.5)
Geo_ALPACA = ROOT.TH2D("Geo_ALPACA", "", 45, -2.5,0.4,30,-7.1,-5.5)

#k= open(pathW+"ALPACA_A.dat","w")
for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        geo=Decimal(float(lm[4])*100.)*Decimal(lm[5])
        if float(lm[4])>0.: Ves_ALPACA.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])), Decimal(float(lm[4])*100.))
        if float(lm[5])>0.: Rec_ALPACA.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])), Decimal(float(lm[5])*100.))
        if geo>0.: Geo_ALPACA.Fill(math.log10(float(lm[0])), math.log10(float(lm[1])), geo)

Ves_ALPACA.SetTitle(";m_{#alpha}(GeV/c^{2});g_{#alpha #gamma}[GeV^{-1}];Vessel Probability(%)")
Ves_ALPACA.SetAxisRange(0.,60.,"Z")
Ves_ALPACA.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Ves_ALPACA.pdf")

Rec_ALPACA.SetTitle(";m_{#alpha}(MeV/c^{2});g_{#alpha #gamma}[GeV^{-1}];Geometric Acceptance(%)")
Rec_ALPACA.SetAxisRange(0.,100.,"Z")
Rec_ALPACA.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Rec_ALPACA.pdf")

Geo_ALPACA.SetTitle(";m_{#alpha}(MeV/c^{2});g_{#alpha #gamma}[GeV^{-1}];Final Acceptance(%)")
Geo_ALPACA.SetAxisRange(0.,60.,"Z")
Geo_ALPACA.Draw("colz")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Geo_ALPACA.pdf")

e_BR,   mu_BR,  tau_BR,         charged_BR,             neutral_BR,             all_BR          = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
e_Vis,  mu_Vis, tau_Vis,        charged_Vis,    neutral_Vis,    all_Vis         = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
e_mass, mu_mass,tau_mass,       charged_mass,   neutral_mass,   all_mass        = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )

m_BR  = ROOT.TMultiGraph()
m_Vis = ROOT.TMultiGraph()
m_BR.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});Branching Ratios")
m_Vis.SetTitle(";m_{#gamma^{D}}(GeV/c^{2});Visible Branching Ratios")


lBR = open(pathW+"BR.dat","r")
l_BR = lBR.readlines()

for l in l_BR:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        if float(lm[1])!=0.: 
                e_BR.append(float(lm[1]))
                e_mass.append(float(lm[0]))
        if float(lm[2])!=0.: 
                mu_BR.append(float(lm[2]))
                mu_mass.append(float(lm[0]))
        if float(lm[3])!=0.: 
                tau_BR.append(float(lm[3]))
                tau_mass.append(float(lm[0]))
        if float(lm[4])!=0.: 
                charged_BR.append(float(lm[4]))
                charged_mass.append(float(lm[0]))
        if float(lm[5])!=0.: 
                neutral_BR.append(float(lm[5]))
                neutral_mass.append(float(lm[0]))
        if float(lm[6])!=0.: 
                all_BR.append(float(lm[6]))
                all_mass.append(float(lm[0]))

eBR                     =ROOT.TGraph(len(e_BR),e_mass,e_BR)
muBR            =ROOT.TGraph(len(mu_BR),mu_mass,mu_BR)
tauBR           =ROOT.TGraph(len(tau_BR),tau_mass,tau_BR)
chargedBR       =ROOT.TGraph(len(charged_BR),charged_mass,charged_BR)
neutralBR       =ROOT.TGraph(len(neutral_BR),neutral_mass,neutral_BR)
allBR           =ROOT.TGraph(len(all_BR),all_mass,all_BR)

eBR.SetTitle("#gamma^{D} #rightarrow e^{ + } + e^{ - }")
muBR.SetTitle("#gamma^{D} #rightarrow #mu^{ + } + #mu^{-}")
tauBR.SetTitle("#gamma^{D} #rightarrow #tau^{ + } + #tau^{ - }")
chargedBR.SetTitle("#gamma^{D} #rightarrow n #times (h^{ + } + h^{ - })")
neutralBR.SetTitle("#gamma^{D} #rightarrow n #times h^{0}")
allBR.SetTitle("All")

allBR .SetMarkerColor(1)
eBR.SetMarkerColor(2)
tauBR.SetMarkerColor(8)
muBR.SetMarkerColor(4)
chargedBR.SetMarkerColor(6)
neutralBR.SetMarkerColor(13) 

eBR .SetMarkerStyle(24)
muBR.SetMarkerStyle(40)
tauBR.SetMarkerStyle(28)
chargedBR.SetMarkerStyle(30)
neutralBR.SetMarkerStyle(32)
allBR.SetMarkerStyle(46) 

allBR .SetLineColor(1)
eBR.SetLineColor(2)
tauBR.SetLineColor(8)
muBR.SetLineColor(4)
chargedBR.SetLineColor(6)
neutralBR.SetLineColor(13)

eBR .SetDrawOption("ALP")
muBR.SetDrawOption("ALP")
tauBR.SetDrawOption("ALP")
chargedBR.SetDrawOption("ALP")
neutralBR.SetDrawOption("ALP")
allBR.SetDrawOption("ALP")

m_BR.Add(eBR)
m_BR.Add(muBR)
m_BR.Add(tauBR)
m_BR.Add(chargedBR)
m_BR.Add(neutralBR)
m_BR.Add(allBR)

m_BR.Draw("ALP")
xaxis=m_BR.GetXaxis()
xaxis.SetLimits(0.002,10.1)
c.SetLogx()
#c.BuildLegend()
Leg=c.BuildLegend(0.12,0.65,0.32,.85)
Leg.Draw("L")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"BR.pdf")
lBR.close()

lvis = open(pathW+"BR_vis.dat","r")
l_Vis = lvis.readlines()

for l in l_Vis:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        if float(lm[1])!=0.: 
                e_Vis.append(float(lm[1]))
                e_mass.append(float(lm[0]))
        if float(lm[2])!=0.: 
                mu_Vis.append(float(lm[2]))
                mu_mass.append(float(lm[0]))
        if float(lm[3])!=0.: 
                tau_Vis.append(float(lm[3]))
                tau_mass.append(float(lm[0]))
        if float(lm[4])!=0.: 
                charged_Vis.append(float(lm[4]))
                charged_mass.append(float(lm[0]))
        if float(lm[6])!=0.: 
                all_Vis.append(float(lm[6]))
                all_mass.append(float(lm[0]))

eVis        =ROOT.TGraph(len(e_Vis),e_mass,e_Vis)
muVis       =ROOT.TGraph(len(mu_Vis),mu_mass,mu_Vis)
tauVis      =ROOT.TGraph(len(tau_Vis),tau_mass,tau_Vis)
chargedVis  =ROOT.TGraph(len(charged_Vis),charged_mass,charged_Vis)
allVis      =ROOT.TGraph(len(all_Vis),all_mass,all_Vis)

eVis.SetTitle("#gamma^{D} #rightarrow e^{ + } + e^{ - }")
muVis.SetTitle("#gamma^{D} #rightarrow #mu^{ + } + #mu^{-}")
tauVis.SetTitle("#gamma^{D} #rightarrow #tau^{ + } + #tau^{ - }")
chargedVis.SetTitle("#gamma^{D} #rightarrow n #times (h^{ + } + h^{ - })")
allVis.SetTitle("All")

allVis .SetMarkerColor(1)
eVis.SetMarkerColor(2)
tauVis.SetMarkerColor(8)
muVis.SetMarkerColor(4)
chargedVis.SetMarkerColor(6)

eVis .SetMarkerStyle(24)
muVis.SetMarkerStyle(40)
tauVis.SetMarkerStyle(28)
chargedVis.SetMarkerStyle(30)
allVis.SetMarkerStyle(46) 

allVis .SetLineColor(1)
eVis.SetLineColor(2)
tauVis.SetLineColor(8)
muVis.SetLineColor(4)
chargedVis.SetLineColor(6)

eVis .SetDrawOption("ALP")
muVis.SetDrawOption("ALP")
tauVis.SetDrawOption("ALP")
chargedVis.SetDrawOption("ALP")
allVis.SetDrawOption("ALP")

m_Vis.Add(eVis)
m_Vis.Add(muVis)
m_Vis.Add(tauVis)
m_Vis.Add(chargedVis)
m_Vis.Add(allVis)

m_Vis.Draw("ALP")
xaxis=m_Vis.GetXaxis()
xaxis.SetLimits(0.002,10.1)
c.SetLogx()
#c.BuildLegend()
Leg=c.BuildLegend(0.12,0.65,0.32,.85)
Leg.Draw("L")
ROOT.gStyle.SetOptStat(0000)
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"Vis.pdf")
lvis.close()

nt.Write()
nt.Print()
nt.Close()
