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

#c1.SetLogx()
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
        dataC = pd.read_csv(pathR+'2100,1.2554921453200000E-13
2200,1.0177145885320000E-13
2300,9.0147564922700000E-14
2400,7.6507145740520000E-14
2500,6.3341741562690000E-14
2600,5.6101955581679900E-14
2700,5.0889712811019900E-14
2800,4.8415675996260000E-14
2900,4.1224367400800000E-14
3000,3.5203836881309900E-14
3100,3.3305349851200000E-14
3200,2.8915053487180000E-14
3300,2.6192321928980000E-14
3400,2.3385679855709900E-14
3500,2.1023353845090000E-14
3600,1.6732699049360000E-14
3700,1.4485183457710000E-14
3800,1.1129765141000000E-14
3900,8.7421448880820000E-15
4000,6.6849782716160000E-15
4100,4.1845619334390000E-15
4100,2.4213858490510000E-15
4000,1.6922635492560000E-15
3900,1.6460082211189900E-15
3800,1.4554322702820000E-15
3700,1.3719617389360000E-15
3600,1.3380489451690000E-15
3500,1.2971852447289900E-15
3400,1.2344855194180000E-15
3300,1.1531186599020000E-15
3200,1.1250349129050000E-15
3100,1.0915710272540000E-15
3000,1.1058494340160000E-15
2900,9.9529999067390000E-16
2800,9.5376196665380000E-16
2700,9.0017303634050000E-16
2600,8.4115494033309900E-16
2500,7.9496132376990000E-16
2400,7.3230500100530000E-16
2300,7.1686932958080000E-16
2200,6.7495061904460000E-16comb_Rate1.csv', header=None)
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
        #data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate2.csv', header=None)
        #dataC = pd.read_csv(pathR+'comb_Rate2.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate2.csv', header=None)
    else:
        data1 = pd.read_csv(pathR+'qcd_Rate1.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate1.csv', header=None)
        #data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate1.csv', header=None)
        #dataC = pd.read_csv(pathR+'comb_Rate1.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate1.csv', header=None)
    d20, d21 = array('f'), array('f') 
    #d30, d31 = array('f'), array('f')
    d10, d11 = array('f'), array('f')
    #dC0, dC1 = array('f'), array('f')
    d00, d01 = array('f'), array('f')
    d40, d41 = array('f'), array('f')
    for i in data1[0]: d10.append(float(i))
    for i in data1[1]: d11.append(float(i))
    for i in data2[0]: d20.append(float(i))
    for i in data2[1]: d21.append(float(i))
    #for i in data3[0]: d30.append(float(i))
    #for i in data3[1]: d31.append(float(i))
    for i in data4[0]: d40.append(float(i))
    for i in data4[1]: d41.append(float(i))
    for i in dataO[0]: d00.append(float(i))
    for i in dataO[1]: d01.append(float(i))
    #for i in dataC[0]: dC0.append(float(i))
    #for i in dataC[1]: dC1.append(float(i))
    p2=ROOT.TGraph(len(data2[1]),d20, d21)
    p4=ROOT.TGraph(len(data4[1]),d40, d41)
    #p3=ROOT.TGraph(len(data3[1]),d30, d31)
    p1=ROOT.TGraph(len(data1[1]),d10, d11)
    #pC=ROOT.TGraph(len(dataC[1]),dC0, dC1)
    p0=ROOT.TGraph(len(dataO[1]),d00, d01)
    p4.SetLineColor(30)
    #p3.SetLineColor(3)
    p2.SetLineColor(2)
    p1.SetLineColor(4)
    #pC.SetLineColor(1)
    p0.SetLineColor(9)
    #pC.SetLineWidth(3)
    p0.SetLineWidth(3)
    p4.SetLineWidth(3)
    #p3.SetLineWidth(3)
    p2.SetLineWidth(3)
    p1.SetLineWidth(3)
    p4.SetTitle("pbrem old") 
    #p3.SetTitle("pbrem new")
    p2.SetTitle("meson")
    p1.SetTitle("QCD")
    #pC.SetTitle("Combined")
    p0.SetTitle("Combined-old_Pbrem")
    mulgr.Add(p1)
    mulgr.Add(p2)
    #mulgr.Add(p3)
    mulgr.Add(p4)
    mulgr.Add(pC)
    #mulgr.Add(p0)

elif comp=="meson":
    if leptophilic:
        data1 = pd.read_csv(pathR+'meson1_Rate2.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate2.csv', header=None)
    else:
        #data1 = pd.read_csv(pathR+'meson1_Rate1.csv', header=None)
        data1 = pd.read_csv(pathR+'meson.csv', header=None)
        data2 = pd.read_csv(pathR+'meson2_Rate1.csv', header=None)
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
