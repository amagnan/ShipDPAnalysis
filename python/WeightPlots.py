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
nt=ROOT.TFile(pathW+"WeightPlots.root","RECREATE")
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

pbrem_weight = open(pathW+"pbrem_weight.dat","r")
l_pbrem_weight = pbrem_weight.readlines()
l1=l_pbrem_weight
mc, ves, rec= array( 'd' ),  array( 'd' ),  array( 'd' )

for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        mc.append(float(lm[5]))
        ves.append(float(lm[6]))
        rec.append(float(lm[7]))

Tru_pbrem       = ROOT.TH2F("Tru_pbrem", "", 100, min(mc),max(mc),100,-19.,-6.)
Ves_pbrem       = ROOT.TH2F("Ves_pbrem", "", 100, min(mc),max(mc),100,-19.,-6.)
Rec_pbrem       = ROOT.TH2F("Rec_pbrem", "", 100, min(mc),max(mc),100,-19.,-6.)

for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        Tru_pbrem.Fill(float(lm[5]),math.log10(float(lm[1])**2.))
        Ves_pbrem.Fill(float(lm[6]),math.log10(float(lm[1])**2.))
        Rec_pbrem.Fill(float(lm[7]),math.log10(float(lm[1])**2.))

Ves_pbrem.SetTitle(";Vessel weights;#varepsilon^{2};")
Ves_pbrem.Draw("colz")
Ves_pbrem.SetAxisRange(0.,25.,"Z")
c.SetLogx()
ROOT.gStyle.SetOptStat(2222)
box =Ves_pbrem.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"VesW_pbrem.pdf")

Rec_pbrem.SetTitle(";Final weights;#varepsilon^{2};")
Rec_pbrem.Draw("colz")
Rec_pbrem.SetAxisRange(0.,25.,"Z")
c.SetLogx()
ROOT.gStyle.SetOptStat(2222)
box =Rec_pbrem.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"RecW_pbrem.pdf")

Tru_pbrem.SetTitle(";Total weights;#varepsilon^{2};")
Tru_pbrem.Draw("colz")
Tru_pbrem.SetAxisRange(0.,25.,"Z")
c.SetLogx()
ROOT.gStyle.SetOptStat(2222)
box =Tru_pbrem.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"TruW_pbrem.pdf")

meson_weight = open(pathW+"meson_weight.dat","r")
l_meson_weight = meson_weight.readlines()
l1=l_meson_weight
mc, ves, rec= array( 'd' ),  array( 'd' ),  array( 'd' )

for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        mc.append(float(lm[5]))
        ves.append(float(lm[6]))
        rec.append(float(lm[7]))

Tru_meson       = ROOT.TH2F("Tru_meson", "", 100, min(mc),max(mc),100,-19.,-6.)
Ves_meson       = ROOT.TH2F("Ves_meson", "", 100, min(mc),max(mc),100,-19.,-6.)
Rec_meson       = ROOT.TH2F("Rec_meson", "", 100, min(mc),max(mc),100,-19.,-6.)

for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        Tru_meson.Fill(float(lm[5]),math.log10(float(lm[1])**2.))
        Ves_meson.Fill(float(lm[6]),math.log10(float(lm[1])**2.))
        Rec_meson.Fill(float(lm[7]),math.log10(float(lm[1])**2.))

Ves_meson.SetTitle(";Vessel weights;#varepsilon^{2};")
Ves_meson.Draw("colz")
Ves_meson.SetAxisRange(0.,25.,"Z")
c.SetLogx()
ROOT.gStyle.SetOptStat(2222)
box = Ves_meson.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"VesW_meson.pdf")

Rec_meson.SetTitle(";Final weights;#varepsilon^{2};")
Rec_meson.Draw("colz")
Rec_meson.SetAxisRange(0.,25.,"Z")
c.SetLogx()
ROOT.gStyle.SetOptStat(2222)
box = Rec_meson.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"RecW_meson.pdf")

Tru_meson.SetTitle(";Total weights;#varepsilon^{2};")
Tru_meson.Draw("colz")
Tru_meson.SetAxisRange(0.,25.,"Z")
c.SetLogx()
ROOT.gStyle.SetOptStat(2222)
box = Tru_meson.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"TruW_meson.pdf")

qcd_weight = open(pathW+"qcd_weight.dat","r")
l_qcd_weight = qcd_weight.readlines()
l1=l_qcd_weight

mc, ves, rec= array( 'd' ),  array( 'd' ),  array( 'd' )

for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        mc.append(float(lm[5]))
        ves.append(float(lm[6]))
        rec.append(float(lm[7]))

Tru_qcd     = ROOT.TH2F("Tru_qcd", "", 100, min(mc),max(mc),100,-18.,-11.)
Ves_qcd     = ROOT.TH2F("Ves_qcd", "", 100, min(mc),max(mc),100,-18.,-11.)
Rec_qcd     = ROOT.TH2F("Rec_qcd", "", 100, min(mc),max(mc),100,-18.,-11.)

for l in l1:
        lm = l.replace('\n','')
        lm = lm.split(' ')
        Tru_qcd.Fill(float(lm[5]),math.log10(float(lm[1])**2.))
        Ves_qcd.Fill(float(lm[6]),math.log10(float(lm[1])**2.))
        Rec_qcd.Fill(float(lm[7]),math.log10(float(lm[1])**2.))

Ves_qcd.SetTitle(";Vessel weights;#varepsilon^{2};")
Ves_qcd.Draw("colz")
Ves_qcd.SetAxisRange(0.,25.,"Z")
ROOT.gStyle.SetOptStat(2222)
c.SetLogx()
box = Ves_qcd.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"VesW_qcd.pdf")

Rec_qcd.SetTitle(";Final weights;#varepsilon^{2};")
Rec_qcd.Draw("colz")
Rec_qcd.SetAxisRange(0.,25.,"Z")
ROOT.gStyle.SetOptStat(2222)
c.SetLogx()
box =Rec_qcd.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"RecW_qcd.pdf")

Tru_qcd.SetTitle(";Total weights;#varepsilon^{2};")
Tru_qcd.Draw("colz")
Tru_qcd.SetAxisRange(0.,25.,"Z")
ROOT.gStyle.SetOptStat(2222)
c.SetLogx()
box =Tru_qcd.GetListOfFunctions().FindObject("stats")
box.SetY2NDC(0.30)
box.SetY1NDC(0.15)
box.SetX2NDC(0.65)
box.SetX1NDC(0.80)
box.Draw()
c.Modified()
c.Update()
c.Write()
c.Print(pathW+"TruW_qcd.pdf")

nt.Write()
nt.Print()
nt.Close()
