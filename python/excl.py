import pandas as pd
import numpy as np
import os,sys,getopt
from array import array

leptophilic = 0
comp = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], "l:c:", ["leptophilic=","comp="])

except getopt.GetoptError:
    print 'exit the system'
    sys.exit()

for o,a in opts:
    if o in ('-l',): leptophilic = a
    if o in ('-c',): comp = a

pathR = "../Exclusion/"

import ROOT
#For the axis titles:
ROOT.gStyle.SetTitleColor(1, "XYZ")
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetTitleSize(0.055, "XYZ")
ROOT.gStyle.SetTitleXOffset(0.9)
ROOT.gStyle.SetTitleYOffset(0.95)

#For the axis labels:
ROOT.gStyle.SetLabelColor(1, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
ROOT.gStyle.SetLabelSize(0.05, "XYZ")

#For the axis:
ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetStripDecimals(True)
ROOT.gStyle.SetTickLength(0.02, "XYZ")
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetPadTickX(1)  #To get tick marks on the opposite side of the frame
ROOT.gStyle.SetPadTickY(1)

c1  = ROOT.TCanvas('c1', '',1920,1080)

c1.SetLogx()
c1.SetLogy()


mulgr = ROOT.TMultiGraph()

if comp==0:
    if leptophilic:
        data1 = pd.read_csv(pathR+'qcd_Rate2.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate2.csv', header=None)
        data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate2.csv', header=None)
    else:
        data1 = pd.read_csv(pathR+'qcd_Rate1.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate1.csv', header=None)
        data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate1.csv', header=None)
    d20, d21 = array('f'), array('f') 
    d30, d31 = array('f'), array('f')
    d10, d11 = array('f'), array('f')
    dC0, dC1 = array('f'), array('f')
    for i in data1[0]: d10.append(float(i))
    for i in data1[1]: d11.append(float(i))
    for i in data2[0]: d20.append(float(i))
    for i in data2[1]: d21.append(float(i))
    for i in data3[0]: d30.append(float(i))
    for i in data3[1]: d31.append(float(i))
    for i in dataC[0]: dC0.append(float(i))
    for i in dataC[1]: dC1.append(float(i))
    p2=ROOT.TGraph(len(d21),d20,d21) 
    p3=ROOT.TGraph(len(data3[1]),d30,d31)
    p1=ROOT.TGraph(len(data1[1]),d10,d11)
    pC=ROOT.TGraph(len(dataC[1]),dC0,dC1)
    p3.SetLineWidth(3)
    p2.SetLineWidth(3)
    p1.SetLineWidth(3)
    pC.SetLineWidth(3)
    p3.SetLineColor(3)
    p2.SetLineColor(2)
    p1.SetLineColor(4)
    pC.SetLineColor(1)
    p3.SetTitle("pbrem")
    p2.SetTitle("meson")
    p1.SetTitle("QCD")
    pC.SetTitle("Combined")
    mulgr.Add(p1)
    mulgr.Add(p2)
    mulgr.Add(p3)
    mulgr.Add(pC) 

elif comp=="pbrem":
    if leptophilic:
        data1 = pd.read_csv(pathR+'qcd_Rate2.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate2.csv', header=None)
        data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate2.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate2.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate2.csv', header=None)
    else:
        data1 = pd.read_csv(pathR+'qcd_Rate1.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate1.csv', header=None)
        data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate1.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate1.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate1.csv', header=None)
    d20, d21 = array('f'), array('f') 
    d30, d31 = array('f'), array('f')
    d10, d11 = array('f'), array('f')
    dC0, dC1 = array('f'), array('f')
    d00, d01 = array('f'), array('f')
    d40, d41 = array('f'), array('f')
    for i in data1[0]: d10.append(float(i))
    for i in data1[1]: d11.append(float(i))
    for i in data2[0]: d20.append(float(i))
    for i in data2[1]: d21.append(float(i))
    for i in data3[0]: d30.append(float(i))
    for i in data3[1]: d31.append(float(i))
    for i in data4[0]: d40.append(float(i))
    for i in data4[1]: d41.append(float(i))
    for i in dataO[0]: d00.append(float(i))
    for i in dataO[1]: d01.append(float(i))
    for i in dataC[0]: dC0.append(float(i))
    for i in dataC[1]: dC1.append(float(i))
    p2=ROOT.TGraph(len(data2[1]),d20, d21)
    p4=ROOT.TGraph(len(data4[1]),d40, d41)
    p3=ROOT.TGraph(len(data3[1]),d30, d31)
    p1=ROOT.TGraph(len(data1[1]),d10, d11)
    pC=ROOT.TGraph(len(dataC[1]),dC0, dC1)
    p0=ROOT.TGraph(len(dataO[1]),d00, d01)
    p4.SetLineColor(30)
    p3.SetLineColor(3)
    p2.SetLineColor(2)
    p1.SetLineColor(4)
    pC.SetLineColor(1)
    p0.SetLineColor(9)
    pC.SetLineWidth(3)
    p0.SetLineWidth(3)
    p4.SetLineWidth(3)
    p3.SetLineWidth(3)
    p2.SetLineWidth(3)
    p1.SetLineWidth(3)
    p4.SetTitle("pbrem old") 
    p3.SetTitle("pbrem new")
    p2.SetTitle("meson")
    p1.SetTitle("QCD")
    pC.SetTitle("Combined")
    p0.SetTitle("Combined-old_Pbrem")
    mulgr.Add(p1)
    mulgr.Add(p2)
    mulgr.Add(p3)
    mulgr.Add(p4)
    mulgr.Add(pC)
    mulgr.Add(p0)

elif comp=="meson":
    if leptophilic:
        data1 = pd.read_csv(pathR+'meson1_Rate2.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate2.csv', header=None)
    else:
        data1 = pd.read_csv(pathR+'meson1_Rate1.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate1.csv', header=None)
    d20, d21 = array('f'), array('f') 
    d10, d11 = array('f'), array('f')
    for i in data1[0]: d10.append(float(i))
    for i in data1[1]: d11.append(float(i))
    for i in data2[0]: d20.append(float(i))
    for i in data2[1]: d21.append(float(i))
    p2=ROOT.TGraph(len(data2[1]), d20, d21)
    p1=ROOT.TGraph(len(data1[1]), d10, d11)
    p2.SetLineColor(2)
    p1.SetLineColor(40)
    p2.SetLineWidth(3)
    p1.SetLineWidth(3)
    p2.SetTitle("meson new")
    p1.SetTitle("meson old")
    p2.SetDrawOption("AL")
    p1.SetDrawOption("AL")
    mulgr.Add(p1)
    mulgr.Add(p2)

elif comp=="pbremOnly":
    if leptophilic:
        data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate2.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate2.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate2.csv', header=None)
    else:
        data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate1.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate1.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate1.csv', header=None)
    d40, d41 = array('f'), array('f') 
    d30, d31 = array('f'), array('f')
    for i in data3[0]: d30.append(float(i))
    for i in data3[1]: d31.append(float(i))
    for i in data4[0]: d40.append(float(i))
    for i in data4[1]: d41.append(float(i))
    p4=ROOT.TGraph(len(data4[1]), d40, d41)
    p3=ROOT.TGraph(len(data3[1]), d30, d31)
    p4.SetLineColor(30)
    p3.SetLineColor(3)
    p4.SetLineWidth(3)
    p3.SetLineWidth(3)
    p4.SetTitle("pbrem old") 
    p3.SetTitle("pbrem new")
    mulgr.Add(p4)
    mulgr.Add(p3)

mulgr.SetTitle(";Mass(MeV);#epsilon^{ 2 }")
mulgr.Draw("AL")
c1.BuildLegend(0.65,0.65,0.9,0.9)
c1.Modify()
c1.Update()
c1.SaveAs(pathR+"Excl_"+str(leptophilic)+"_"+str(comp)+".pdf")
c1.Print(pathR+"Excl_"+str(leptophilic)+"_"+str(comp)+".root")
