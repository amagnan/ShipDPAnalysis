import ROOT as r
import os,sys,getopt
import rootUtils as ut
import shipunit as u
import darkphoton
from ShipGeoConfig import ConfigRegistry 
from decorators import *
from rootpyPickler import Unpickler
from array import array
import shipRoot_conf
import dpProductionRates as dputil
import math as m
import numpy as np

shipRoot_conf.configure()

PDG = r.TDatabasePDG.Instance()

try:
    opts, args = getopt.getopt(sys.argv[1:], "f:g:", ["nEvents=","geoFile="])

except getopt.GetoptError:
    print 'no file'
    sys.exit()

for o,a in opts:
    if o in ('-f',): inputFile = a
    if o in ('-g', '--geoFile',): geoFile = a

tmp=inputFile.replace("/eos/experiment/ship/data/DarkPhoton/PBC-June-3/sim/","")
tmp1=tmp.replace('.root','')
tmp2=tmp1.replace('mass','')
tmp3=tmp2.replace('eps','')     
out=tmp3.split('_') 
pro=out[0]
mass_mc=float(out[1])
eps=float(out[2])


eosship =  r.gSystem.Getenv("EOSSHIP")
eospath = eosship+inputFile
f = r.TFile.Open(eospath)
sTree=f.cbmsim
eospath = eosship+geoFile
fgeo = r.TFile.Open(eospath)

sGeo = r.gGeoManager
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')

dy = ShipGeo.Yheight/u.m

import shipDet_conf

run = r.FairRunSim()
run.SetName('TGeant4')
run.SetOutputFile('dummy') 
run.SetUserConfig('g4Config_basic.C')
rtdb = run.GetRuntimeDb()

modules = shipDet_conf.configure(run,ShipGeo)
run.Init()

import geomGeant4

if hasattr(ShipGeo.Bfield,"fieldMap"):
    fieldMaker = geomGeant4.addVMCFields(ShipGeo.Bfield.fieldMap, ShipGeo.Bfield.z, True)

hcalGeoFile = ShipGeo.hcal.File
ecalGeoFile = ShipGeo.ecal.File
dy = ShipGeo.Yheight/u.m
MeasCut=25

geoMat =  r.genfit.TGeoMaterialInterface()
r.genfit.MaterialEffects.getInstance().init(geoMat)
bfield = r.genfit.FairShipFields()
fM = r.genfit.FieldManager.getInstance()
fM.init(bfield)


h={}

ut.bookHist(h,'DP','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)

ut.bookHist(h,'DP_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_tau','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_tau','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_mix','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_mix','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_single','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_single','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)

ut.bookHist(h,'DPang_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_tau','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_tau','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_mix','invariant Mass(GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_mix','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_single','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_single','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)

ut.bookHist(h,'DPves_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_tau','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_tau','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_mix','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_mix','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_single','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_single','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)

ut.bookHist(h,'eff_P','',50,0.,250.)
ut.bookHist(h,'eff_Pt','',50,0.,8.)
ut.bookHist(h,'eff_Eta','',50,2.,8.)
ut.bookHist(h,'eff_Pw','',50,0.,250.)
ut.bookHist(h,'eff_Ptw','',50,0.,8.)
ut.bookHist(h,'eff_Etaw','',50,2.,8.)

ut.bookHist(h,'eff_P_e','',50,0.,250.)
ut.bookHist(h,'eff_Pt_e','',50,0.,8.)
ut.bookHist(h,'eff_Eta_e','',50,2.,8.)
ut.bookHist(h,'eff_P_mu','',50,0.,250.)
ut.bookHist(h,'eff_Pt_mu','',50,0.,8.)
ut.bookHist(h,'eff_Eta_mu','',50,2.,8.)
ut.bookHist(h,'eff_P_tau','',50,0.,250.)
ut.bookHist(h,'eff_Pt_tau','',50,0.,8.)
ut.bookHist(h,'eff_Eta_tau','',50,2.,8.)
ut.bookHist(h,'eff_P_pi','',50,0.,250.)
ut.bookHist(h,'eff_Pt_pi','',50,0.,8.)
ut.bookHist(h,'eff_Eta_pi','',50,2.,8.)
ut.bookHist(h,'eff_P_ka','',50,0.,250.)
ut.bookHist(h,'eff_Pt_ka','',50,0.,8.)
ut.bookHist(h,'eff_Eta_ka','',50,2.,8.)
ut.bookHist(h,'eff_P_oth','',50,0.,250.)
ut.bookHist(h,'eff_Pt_oth','',50,0.,8.)
ut.bookHist(h,'eff_Eta_oth','',50,2.,8.)


ut.bookHist(h,'eff_P0w','',50,0.,250.)
ut.bookHist(h,'eff_P1w','',50,0.,250.)
ut.bookHist(h,'eff_P2w','',50,0.,250.)
ut.bookHist(h,'eff_Pt0w','',50,0.,8.)
ut.bookHist(h,'eff_Pt1w','',50,0.,8.)
ut.bookHist(h,'eff_Pt2w','',50,0.,8.)
ut.bookHist(h,'eff_Eta0w','',50,2.,8.)
ut.bookHist(h,'eff_Eta1w','',50,2.,8.)
ut.bookHist(h,'eff_Eta2w','',50,2.,8.)

ut.bookHist(h,'eff_P0','',50,0.,250.)
ut.bookHist(h,'eff_P1','',50,0.,250.)
ut.bookHist(h,'eff_P2','',50,0.,250.)
ut.bookHist(h,'eff_Pt0','',50,0.,8.)
ut.bookHist(h,'eff_Pt1','',50,0.,8.)
ut.bookHist(h,'eff_Pt2','',50,0.,8.)
ut.bookHist(h,'eff_Eta0','',50,2.,8.)
ut.bookHist(h,'eff_Eta1','',50,2.,8.)
ut.bookHist(h,'eff_Eta2','',50,2.,8.)

ut.bookHist(h,'eff_P0_e','',50,0.,250.)
ut.bookHist(h,'eff_P1_e','',50,0.,250.)
ut.bookHist(h,'eff_P2_e','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_e','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_e','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_e','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_e','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_e','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_e','',50,2.,8.)

ut.bookHist(h,'eff_P0_mu','',50,0.,250.)
ut.bookHist(h,'eff_P1_mu','',50,0.,250.)
ut.bookHist(h,'eff_P2_mu','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_mu','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_mu','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_mu','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_mu','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_mu','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_mu','',50,2.,8.)

ut.bookHist(h,'eff_P0_tau','',50,0.,250.)
ut.bookHist(h,'eff_P1_tau','',50,0.,250.)
ut.bookHist(h,'eff_P2_tau','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_tau','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_tau','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_tau','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_tau','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_tau','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_tau','',50,2.,8.)

ut.bookHist(h,'eff_P0_pi','',50,0.,250.)
ut.bookHist(h,'eff_P1_pi','',50,0.,250.)
ut.bookHist(h,'eff_P2_pi','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_pi','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_pi','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_pi','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_pi','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_pi','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_pi','',50,2.,8.)

ut.bookHist(h,'eff_P0_ka','',50,0.,250.)
ut.bookHist(h,'eff_P1_ka','',50,0.,250.)
ut.bookHist(h,'eff_P2_ka','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_ka','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_ka','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_ka','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_ka','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_ka','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_ka','',50,2.,8.)

ut.bookHist(h,'eff_P0_oth','',50,0.,250.)
ut.bookHist(h,'eff_P1_oth','',50,0.,250.)
ut.bookHist(h,'eff_P2_oth','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_oth','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_oth','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_oth','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_oth','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_oth','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_oth','',50,2.,8.)

ut.bookHist(h,'eff_P0_mix','',50,0.,250.)
ut.bookHist(h,'eff_P1_mix','',50,0.,250.)
ut.bookHist(h,'eff_P2_mix','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_mix','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_mix','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_mix','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_mix','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_mix','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_mix','',50,2.,8.)

h['eff_P1_mu'].Sumw2()
h['eff_Pt1_mu'].Sumw2()
h['eff_Eta1_mu'].Sumw2()
h['eff_P2_mu'].Divide(h['eff_P1_mu']*2e20, h['eff_P0_mu'],1.,1.,"B")
h['eff_Pt2_mu'].Divide(h['eff_Pt1_mu']*2e20, h['eff_Pt0_mu'],1.,1.,"B")
h['eff_Eta2_mu'].Divide(h['eff_Eta1_mu']*2e20,h['eff_Eta0_mu'],1.,1.,"B")
h['eff_P2_mu'].Draw("E")
h['eff_P2_mu'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_mu'].Draw("E")
h['eff_Pt2_mu'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_mu'].Draw("E")
h['eff_Eta2_mu'].SetTitle("PseudoRapidity ;#eta;Eff.")

def dist2InnerWall(X,Y,Z):
    dist = 0
    node = sGeo.FindNode(X,Y,Z)
    if ShipGeo.tankDesign < 5:
        if not 'cave' in node.GetName(): return dist  # TP 
    else:
        if not 'decayVol' in node.GetName(): return dist
    start = array('d',[X,Y,Z])
    nsteps = 8
    dalpha = 2*r.TMath.Pi()/nsteps
    rsq = X**2+Y**2
    minDistance = 100 *u.m
    for n in range(nsteps):
        alpha = n * dalpha
        sdir  = array('d',[r.TMath.Sin(alpha),r.TMath.Cos(alpha),0.])
        node = sGeo.InitTrack(start, sdir)
        nxt = sGeo.FindNextBoundary()
        if ShipGeo.tankDesign < 5 and nxt.GetName().find('I')<0: return 0    
        distance = sGeo.GetStep()
        if distance < minDistance  : minDistance = distance
    return minDistance

def isInFiducial(X,Y,Z):
    if Z > ShipGeo.TrackStation1.z : return False
    if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
    if dist2InnerWall(X,Y,Z)<5*u.cm: return False
    return True 

def checkFiducialVolume(vtx,a):
    x, y, z = vtx.X(), vtx.Y(), vtx.Z()
    x_dprime = 1000.
    y_dprime = 768.2500
    z_dprime0 = 1349.7500
    x_prime0=x+z_dprime0*a[0]
    y_prime0=y+z_dprime0*a[1]
    if x_prime0<=abs(x_dprime) and y_prime0<=abs(y_dprime): return True

def ang_muon(vtx,a):
    x, y, z = vtx.X(), vtx.Y(), vtx.Z()
    x_dprime = 300.
    y_dprime = 600.
    z_dprime0 = 3918.6000
    x_prime0=x+z_dprime0*a[0]
    y_prime0=y+z_dprime0*a[1]
    if x_prime0<=abs(x_dprime) and y_prime0<=abs(y_dprime): return True

def ang_ecal(vtx,a):
    x, y, z = vtx.X(), vtx.Y(), vtx.Z()
    x_dprime = 265.
    y_dprime = 530.
    z_dprime0 = 3621.5501
    x_prime0=x+z_dprime0*a[0]
    y_prime0=y+z_dprime0*a[1]
    if x_prime0<=abs(x_dprime) and y_prime0<=abs(y_dprime): return True


def ang_hcal(vtx,a):
    x, y, z = vtx.X(), vtx.Y(), vtx.Z()
    x_dprime = 315.
    y_dprime = 630.
    z_dprime0 = 3666.6001
    x_prime0=x+z_dprime0*a[0]
    y_prime0=y+z_dprime0*a[1]
    if x_prime0<=abs(x_dprime) and y_prime0<=abs(y_dprime): return True



def wgnew(dark):
    flmin=47.505*u.m
    flmax=98.265*u.m
    LS = r.gRandom.Uniform(flmin,flmax)
    p = dark.GetP()
    e = dark.GetEnergy()
    gam = e/r.TMath.Sqrt(e*e-p*p)
    beta = p/e
    DP_instance = darkphoton.DarkPhoton(mass_mc,eps)
    ctau = DP_instance.cTau()
    return r.TMath.Exp(-LS/(beta*gam*ctau))*( (ShipGeo.NuTauTarget.zdim*2.)/(beta*gam*ctau) )

def angle(P,Px,Py,Pz):
    m_x, m_y, m_z   = Px, Py, Pz
    mom=P
    a_x=m.tan(m_x/mom)
    a_y=m.tan(m_y/mom)
    return a_x,a_y

def invariant(S_x,S_y,e_0,e_1,p_0,p_1,T,Ang):
    #def invariant(S_x,S_y,m_0,m_1,p_0,p_1):
    th=r.TMath.Sqrt((S_x[0]-S_x[1])**2.+(S_y[0]-S_y[1])**2.)
    #th=r.TMath.Sqrt(S_x[0]**2.+S_y[0]**2.)+r.TMath.Sqrt(S_x[1]**2.+S_y[1]**2.)
    #print T[0]+T[1],Ang[0]+Ang[1]
    #th=abs(Ang[0]-Ang[1])
    #e_0=m_0**2.+p_0**2.
    #e_1=m_1**2.+p_1**2.
    #return r.TMath.Sqrt(m_0**2.+m_1**2.+2.*(e_0*e_1-p_0*p_1*m.cos(th)))
    return r.TMath.Sqrt((e_0+e_1)**2.-(p_0**2.+p_1**2.+2.*p_0*p_1*m.cos(th)))

def myEventLoop(n,xsw):
    neg,pos=0,0
    dp, e,  mu, tau,    pi, ka, oth =0, 0, 0, 0, 0, 0, 0
    dp_a,   e_a,    mu_a,   tau_a,  pi_a,   ka_a,   oth_a   =0, 0, 0, 0, 0, 0, 0
    dp_v,   e_v,    mu_v,   tau_v,  pi_v,   ka_v,   oth_v   =0, 0, 0, 0, 0, 0, 0
    rc=sTree.GetEntry(n)
    had,vec=[],[]
    mo=-99
    Ang, E,P,M,S_y,S_x,T=[],[],[],[],[],[],[]
    if pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,'meson',sTree.MCTrack[0].GetPdgCode())
    wg = sTree.MCTrack[1].GetWeight()
    if abs(wg)<0.00000001: wg = wgnew(sTree.MCTrack[1])
    for mc,track in enumerate(sTree.MCTrack):
        """if track.GetMotherId()>0 and sTree.MCTrack[track.GetMotherId()].GetPdgCode()!=9900015 and sTree.MCTrack[track.GetMotherId()].GetPdgCode()!=4900023 and (track.GetPdgCode()==9900015 or track.GetPdgCode()==4900023):
            mo=mc
            print mo
        if track.GetMotherId()==mo:"""
        if track.GetMotherId()==1 and (sTree.MCTrack[1].GetPdgCode()==9900015 or sTree.MCTrack[1].GetPdgCode()==4900023):
            if abs(track.GetPdgCode())==11 or abs(track.GetPdgCode())==13 or abs(track.GetPdgCode())==15:
                dp+=1
                h['eff_P'].Fill(track.GetP(),wg)            
                h['eff_Pt'].Fill(track.GetPt(),wg)
                h['eff_Eta'].Fill(track.GetRapidity(),wg)
                X=track.GetStartX()
                Y=track.GetStartY()
                Z=track.GetStartZ()
                Px=track.GetPx()
                Py=track.GetPy()
                Pz=track.GetPz()
                P.append(track.GetP())
                E.append(track.GetEnergy())
                a=angle(track.GetP(),Px,Py,Pz)
                T.append(track.GetStartT())
                Ang.append(m.atan(track.GetPt()/Pz))
                vtx=r.TVector3(X,Y,Z)
                mom=r.TLorentzVector()
                track.Get4Momentum(mom)
                M.append(mom.M())
                #print r.TMath.Sqrt(Px**2.+Py**2.+Pz**2.),r.TMath.Sqrt(mom.X()**2.+mom.Y()**2.+mom.Z()**2.)
                S_x.append(a[0])
                S_y.append(a[1])
                if isInFiducial(X,Y,Z):
                    dp_v+=1
                    h['eff_P0'].Fill(track.GetP())          
                    h['eff_Pt0'].Fill(track.GetPt())
                    h['eff_Eta0'].Fill(track.GetRapidity())
                    if ang_ecal(vtx,a) or ang_muon(vtx,a):
                        dp_a+=1
                        h['eff_P1'].Fill(track.GetP())          
                        h['eff_Pt1'].Fill(track.GetPt())
                        h['eff_Eta1'].Fill(track.GetRapidity())
                if abs(track.GetPdgCode())==11:
                    e+=1
                    h['eff_P_e'].Fill(track.GetP(),wg)          
                    h['eff_Pt_e'].Fill(track.GetPt(),wg)
                    h['eff_Eta_e'].Fill(track.GetRapidity(),wg)
                    if isInFiducial(X,Y,Z):
                        e_v+=1
                        h['eff_P0_e'].Fill(track.GetP())            
                        h['eff_Pt0_e'].Fill(track.GetPt())
                        h['eff_Eta0_e'].Fill(track.GetRapidity())
                        if ang_ecal(vtx,a):
                            e_a+=1
                            h['eff_P1_e'].Fill(track.GetP())            
                            h['eff_Pt1_e'].Fill(track.GetPt())
                            h['eff_Eta1_e'].Fill(track.GetRapidity())
                if abs(track.GetPdgCode())==13:
                    mu+=1
                    h['eff_P_mu'].Fill(track.GetP(),wg)         
                    h['eff_Pt_mu'].Fill(track.GetPt(),wg)
                    h['eff_Eta_mu'].Fill(track.GetRapidity(),wg)
                    if isInFiducial(X,Y,Z):
                        mu_v+=1
                        h['eff_P0_mu'].Fill(track.GetP())           
                        h['eff_Pt0_mu'].Fill(track.GetPt())
                        h['eff_Eta0_mu'].Fill(track.GetRapidity())
                        if ang_muon(vtx,a):
                            mu_a+=1
                            h['eff_P1_mu'].Fill(track.GetP())           
                            h['eff_Pt1_mu'].Fill(track.GetPt())
                            h['eff_Eta1_mu'].Fill(track.GetRapidity())
                if abs(track.GetPdgCode())==15:#dont think that is the right way to do that
                    tau+=1
                    h['eff_P_tau'].Fill(track.GetP(),wg)            
                    h['eff_Pt_tau'].Fill(track.GetPt(),wg)
                    h['eff_Eta_tau'].Fill(track.GetRapidity(),wg)
                    if isInFiducial(X,Y,Z):     
                        tau_v+=1
                        h['eff_P0_tau'].Fill(track.GetP())          
                        h['eff_Pt0_tau'].Fill(track.GetPt())
                        h['eff_Eta0_tau'].Fill(track.GetRapidity())
                        if ang_ecal(vtx,a):
                            tau_a+=1
                            h['eff_P1_tau'].Fill(track.GetP())          
                            h['eff_Pt1_tau'].Fill(track.GetPt())
                            h['eff_Eta1_tau'].Fill(track.GetRapidity())

        if track.GetMotherId()==1 and abs(track.GetPdgCode())<7:
            had.append(mc)
            #print track.GetPdgCode()
            #if abs(track.GetPdgCode())>6: print track.GetPdgCode()#for testing else state
        if dp != 2:
            for daug in had:
                X=track.GetStartX()
                Y=track.GetStartY()
                Z=track.GetStartZ()
                Px=track.GetPx()
                Py=track.GetPy()
                Pz=track.GetPz()
                vtx=r.TVector3(X,Y,Z)
                mom=r.TLorentzVector()
                track.Get4Momentum(mom)
                if track.GetMotherId()==daug:
                    P.append(track.GetP())
                    E.append(track.GetEnergy())
                    M.append(mom.M())
                    Ang.append(m.atan(track.GetPt()/Pz))
                    T.append(track.GetStartT())
                    #print r.TMath.Sqrt(Px**2.+Py**2.+Pz**2.),r.TMath.Sqrt(mom.X()**2.+mom.Y()**2.+mom.Z()**2.)
                    a=angle(track.GetP(),Px,Py,Pz)
                    S_x.append(a[0])
                    S_y.append(a[1])
                    #print n,track.GetPdgCode()
                    if abs(track.GetPdgCode())==211:
                        pi+=1
                        h['eff_P_pi'].Fill(track.GetP(),wg)         
                        h['eff_Pt_pi'].Fill(track.GetPt(),wg)
                        h['eff_Eta_pi'].Fill(track.GetRapidity(),wg)
                        if isInFiducial(X,Y,Z):
                            pi_v+=1
                            h['eff_P0_pi'].Fill(track.GetP())           
                            h['eff_Pt0_pi'].Fill(track.GetPt())
                            h['eff_Eta0_pi'].Fill(track.GetRapidity())
                            if ang_ecal(vtx,a):
                                pi_a+=1
                                h['eff_P1_pi'].Fill(track.GetP())           
                                h['eff_Pt1_pi'].Fill(track.GetPt())
                                h['eff_Eta1_pi'].Fill(track.GetRapidity())
                    elif abs(track.GetPdgCode())==321:
                        ka+=1
                        h['eff_P_ka'].Fill(track.GetP(),wg)         
                        h['eff_Pt_ka'].Fill(track.GetPt(),wg)
                        h['eff_Eta_ka'].Fill(track.GetRapidity(),wg)
                        if isInFiducial(X,Y,Z):
                            ka_v+=1
                            h['eff_P0_ka'].Fill(track.GetP())           
                            h['eff_Pt0_ka'].Fill(track.GetPt())
                            h['eff_Eta0_ka'].Fill(track.GetRapidity())
                            if ang_ecal(vtx,a):
                                ka_a+=1
                                h['eff_P1_ka'].Fill(track.GetP())           
                                h['eff_Pt1_ka'].Fill(track.GetPt())
                                h['eff_Eta1_ka'].Fill(track.GetRapidity())
                    else:
                        #print n,track.GetPdgCode()
                        oth+=1
                        h['eff_P_oth'].Fill(track.GetP(),wg)            
                        h['eff_Pt_oth'].Fill(track.GetPt(),wg)
                        h['eff_Eta_oth'].Fill(track.GetRapidity(),wg)
                        if isInFiducial(X,Y,Z):
                            oth_v+=1
                            h['eff_P0_oth'].Fill(track.GetP())          
                            h['eff_Pt0_oth'].Fill(track.GetPt())
                            h['eff_Eta0_oth'].Fill(track.GetRapidity())
                            if ang_ecal(vtx,a) or ang_hcal(vtx,a):
                                oth_a+=1
                                h['eff_P1_oth'].Fill(track.GetP())          
                                h['eff_Pt1_oth'].Fill(track.GetPt())
                                h['eff_Eta1_oth'].Fill(track.GetRapidity())
                    if abs(track.GetPdgCode())==211 or abs(track.GetPdgCode())==321:
                        if track.GetPdgCode()<0: neg=1
                        if track.GetPdgCode()>0: pos=1
                    dp+=1
                    if abs(track.GetPdgCode())==221 or abs(track.GetPdgCode())==213 or abs(track.GetPdgCode())==223: vec.append(mc)
                    if isInFiducial(X,Y,Z):
                        dp_v+=1
                        h['eff_P0'].Fill(track.GetP())          
                        h['eff_Pt0'].Fill(track.GetPt())
                        h['eff_Eta0'].Fill(track.GetRapidity())
                        if ang_ecal(vtx,a) or ang_hcal(vtx,a) or ang_muon(vtx,a):
                            dp_a+=1
                            h['eff_P1'].Fill(track.GetP())          
                            h['eff_Pt1'].Fill(track.GetPt())
                            h['eff_Eta1'].Fill(track.GetRapidity())


    if dp==2:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1],T,Ang)
        h['eff_Pw'].Fill(sTree.MCTrack[1].GetP(),wg)            
        h['eff_Ptw'].Fill(sTree.MCTrack[1].GetPt(),wg)
        h['eff_Etaw'].Fill(sTree.MCTrack[1].GetRapidity(),wg)
        h['DP'].Fill(mass)
        h['DPW'].Fill(mass,wg*xsw)
        if e==2:  
            h['DP_e'].Fill(mass)
            h['DPW_e'].Fill(mass,wg*xsw)
        if mu==2:
            #print mass              
            h['DP_mu'].Fill(mass)
            h['DPW_mu'].Fill(mass,wg*xsw)
        if tau==2:              
            h['DP_tau'].Fill(mass)
            h['DPW_tau'].Fill(mass,wg*xsw)
        if pi==2:      
            h['DP_pi'].Fill(mass)
            h['DPW_pi'].Fill(mass,wg*xsw)
        if ka==2:
            h['DP_ka'].Fill(mass)       
            h['DPW_ka'].Fill(mass,wg*xsw)
        if oth==2:
            h['DP_oth'].Fill(mass)      
            h['DPW_oth'].Fill(mass,wg*xsw)
        if oth==1:
            h['DP_mix'].Fill(mass)      
            h['DPW_mix'].Fill(mass,wg*xsw)

    if dp==1 and oth==1:
        h['DP_single'].Fill(mass_mc)        
        h['DPW_single'].Fill(mass_mc,wg*xsw)
    if dp_v==1 and oth_v==1:
        h['DPves_single'].Fill(mass_mc)     
        h['DPvesW_single'].Fill(mass_mc,wg*xsw)
    if dp_a==1 and oth_a==1:
        h['DPang_single'].Fill(mass_mc)     
        h['DPangW_single'].Fill(mass_mc,wg*xsw)

    if dp_v==2:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1],T,Ang)
        h['DPves'].Fill(mass)
        h['DPvesW'].Fill(mass,wg*xsw)
        h['eff_P0w'].Fill(sTree.MCTrack[1].GetP(),wg*xsw)           
        h['eff_Pt0w'].Fill(sTree.MCTrack[1].GetPt(),wg*xsw)
        h['eff_Eta0w'].Fill(sTree.MCTrack[1].GetRapidity(),wg*xsw)
        if e_v==2:  
            h['DPves_e'].Fill(mass)
            h['DPvesW_e'].Fill(mass,wg*xsw)
        if mu_v==2:              
            h['DPves_mu'].Fill(mass)
            h['DPvesW_mu'].Fill(mass,wg*xsw)
        if tau_v==2:              
            h['DPves_tau'].Fill(mass)
            h['DPvesW_tau'].Fill(mass,wg*xsw)
        if pi_v==2:             
            h['DPves_pi'].Fill(mass)
            h['DPvesW_pi'].Fill(mass,wg*xsw)
        if ka_v==2:
            h['DPves_ka'].Fill(mass)        
            h['DPvesW_ka'].Fill(mass,wg*xsw)
        if oth_v==2:
            h['DPves_oth'].Fill(mass)       
            h['DPvesW_oth'].Fill(mass,wg*xsw)
        if oth_v==1:
            h['DPves_mix'].Fill(mass)       
            h['DPvesW_mix'].Fill(mass,wg*xsw)

    if dp_a==2:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1],T,Ang)
        h['DPang'].Fill(mass)
        h['DPangW'].Fill(mass,wg*xsw)
        h['eff_P1w'].Fill(sTree.MCTrack[1].GetP(),wg*xsw)           
        h['eff_Pt1w'].Fill(sTree.MCTrack[1].GetPt(),wg*xsw)
        h['eff_Eta1w'].Fill(sTree.MCTrack[1].GetRapidity(),wg*xsw)
        if e_a==2:  
            h['DPang_e'].Fill(mass)
            h['DPangW_e'].Fill(mass,wg*xsw)
        if mu_a==2:              
            h['DPang_mu'].Fill(mass)
            h['DPangW_mu'].Fill(mass,wg*xsw)
        if tau_a==2:              
            h['DPang_tau'].Fill(mass)
            h['DPangW_tau'].Fill(mass,wg*xsw)
        if pi_a==2:              
            h['DPang_pi'].Fill(mass)
            h['DPangW_pi'].Fill(mass,wg*xsw)
        if ka_a==2:
            h['DPang_ka'].Fill(mass)        
            h['DPangW_ka'].Fill(mass,wg*xsw)
        if oth_a==2:
            h['DPang_oth'].Fill(mass)       
            h['DPangW_oth'].Fill(mass,wg*xsw)
        if oth_a==1:
            h['DPang_mix'].Fill(mass)       
            h['DPangW_mix'].Fill(mass,wg*xsw)


nEvents =sTree.GetEntries()
print nEvents
for n in range(nEvents):
    if not pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,pro,0)
    myEventLoop(n,xsw)

o1 = tmp1+"_rec.dat"
o2 = tmp1+"_e.dat"
o3 = tmp1+"_mu.dat"
o4 = tmp1+"_pi.dat"
o5 = tmp1+"_ka.dat"
o6 = tmp1+"_oth.dat"
o7 = tmp1+"_mix.dat"
o8 = tmp1+"_tau.dat"
#o9 = tmp1+"_single.dat"
a=open(o1,'w+')
b=open(o2,'w+')
c=open(o3,'w+')
d=open(o4,'w+')
e=open(o5,'w+')
f=open(o6,'w+')
g=open(o7,'w+')
t=open(o8,'w+')
#s=open(o9,w+)

o11 = tmp1+"_table_rec.dat"
o12 = tmp1+"_table_e.dat"
o13 = tmp1+"_table_mu.dat"
o14 = tmp1+"_table_pi.dat"
o15 = tmp1+"_table_ka.dat"
o16 = tmp1+"_table_oth.dat"
o17 = tmp1+"_table_mix.dat"
o18 = tmp1+"_table_tau.dat"
o19 = tmp1+"_table_single.dat"
a1=open(o11,'w+')
b1=open(o12,'w+')
c1=open(o13,'w+')
d1=open(o14,'w+')
e1=open(o15,'w+')
f1=open(o16,'w+')
g1=open(o17,'w+')
t1=open(o18,'w+')
s1=open(o19,'w+')

if float(h['DPvesW'].Integral())!=0: Acc=float(h['DPangW'].Integral())/float(h['DPvesW'].Integral())
if float(h['DPvesW_e'].Integral())!=0: Acc_e=float(h['DPangW_e'].Integral())/float(h['DPvesW_e'].Integral())
if float(h['DPvesW_mu'].Integral())!=0: Acc_mu=float(h['DPangW_mu'].Integral())/float(h['DPvesW_mu'].Integral())
if float(h['DPvesW_tau'].Integral())!=0: Acc_tau=float(h['DPangW_tau'].Integral())/float(h['DPvesW_tau'].Integral())
if float(h['DPvesW_pi'].Integral())!=0: Acc_pi=float(h['DPangW_pi'].Integral())/float(h['DPvesW_pi'].Integral())
if float(h['DPvesW_ka'].Integral())!=0: Acc_ka=float(h['DPangW_ka'].Integral())/float(h['DPvesW_ka'].Integral())
if float(h['DPvesW_oth'].Integral())!=0: Acc_oth=float(h['DPangW_oth'].Integral())/float(h['DPvesW_oth'].Integral())
if float(h['DPvesW_mix'].Integral())!=0: Acc_mix=float(h['DPangW_mix'].Integral())/float(h['DPvesW_mix'].Integral())

if float(h['DPves'].Integral())!=0: 
    RecW=float(h['DPangW'].Integral())/float(h['DPves'].Integral())*2e20
    a.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW, RecW, Acc)) 
    a.write('\n')
if float(h['DPves_e'].Integral())!=0:
    RecW_e=float(h['DPangW_e'].Integral())/float(h['DPves_e'].Integral())*2e20
    b.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_e, RecW_e, Acc_e))
    b.write('\n')
if float(h['DPves_mu'].Integral())!=0:
    RecW_mu=float(h['DPangW_mu'].Integral())/float(h['DPves_mu'].Integral())*2e20
    c.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_mu, RecW_mu, Acc_mu))
    c.write('\n')
if float(h['DPves_tau'].Integral())!=0:
    RecW_tau=float(h['DPangW_tau'].Integral())/float(h['DPves_tau'].Integral())*2e20
    t.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_tau, RecW_tau, Acc_tau))
    t.write('\n')
if float(h['DPves_pi'].Integral())!=0:
    RecW_pi=float(h['DPangW_pi'].Integral())/float(h['DPves_pi'].Integral())*2e20
    d.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_pi, RecW_pi, Acc_pi))
    d.write('\n')
if float(h['DPves_ka'].Integral())!=0:
    RecW_ka=float(h['DPangW_ka'].Integral())/float(h['DPves_ka'].Integral())*2e20
    e.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_ka, RecW_ka, Acc_ka)) 
    e.write('\n')
if float(h['DPves_oth'].Integral())!=0:
    RecW_oth=float(h['DPangW_oth'].Integral())/float(h['DPves_oth'].Integral())*2e20
    f.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_oth, RecW_oth, Acc_oth))
    f.write('\n')
if float(h['DPves_mix'].Integral())!=0:
    RecW_mix=float(h['DPangW_mix'].Integral())/float(h['DPves_mix'].Integral())*2e20
    g.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_mix, RecW_mix, Acc_mix))
    g.write('\n')


a1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP'].Integral()), float(h['DPves'].Integral()), float(h['DPang'].Integral())))
a1.write('\n')
b1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_e'].Integral()), float(h['DPves_e'].Integral()), float(h['DPang_e'].Integral())))
b1.write('\n')
c1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_mu'].Integral()), float(h['DPves_mu'].Integral()), float(h['DPang_mu'].Integral())))
c1.write('\n')
t1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_tau'].Integral()), float(h['DPves_tau'].Integral()), float(h['DPang_tau'].Integral())))
t1.write('\n')
d1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_pi'].Integral()), float(h['DPves_pi'].Integral()), float(h['DPang_pi'].Integral())))
d1.write('\n')
e1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_ka'].Integral()), float(h['DPves_ka'].Integral()), float(h['DPang_ka'].Integral())))
e1.write('\n')
f1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_oth'].Integral()), float(h['DPves_oth'].Integral()), float(h['DPang_oth'].Integral())))
f1.write('\n')
g1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_mix'].Integral()), float(h['DPves_mix'].Integral()), float(h['DPang_mix'].Integral())))
g1.write('\n')
s1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_mix'].Integral()), float(h['DPves_mix'].Integral()), float(h['DPang_mix'].Integral())))
s1.write('\n')

a.close()
b.close()
c.close()
d.close()
e.close()
f.close()
t.close()
g.close()
#s.close()

a1.close()
b1.close()
c1.close()
d1.close()
e1.close()
f1.close()
t1.close()
g1.close()
s1.close()

h['eff_P1w'].Sumw2()
h['eff_Pt1w'].Sumw2()
h['eff_Eta1w'].Sumw2()
h['eff_P2w'].Divide(h['eff_P1w']*2e20, h['eff_P0w'],1.,1.,"B")
h['eff_Pt2w'].Divide(h['eff_Pt1w']*2e20, h['eff_Pt0w'],1.,1.,"B")
h['eff_Eta2w'].Divide(h['eff_Eta1w']*2e20,h['eff_Eta0w'],1.,1.,"B")
h['eff_P2w'].Draw("E")
h['eff_P2w'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2w'].Draw("E")
h['eff_Pt2w'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2w'].Draw("E")
h['eff_Eta2w'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1'].Sumw2()
h['eff_Pt1'].Sumw2()
h['eff_Eta1'].Sumw2()
h['eff_P2'].Divide(h['eff_P1'], h['eff_P0'],1.,1.,"B")
h['eff_Pt2'].Divide(h['eff_Pt1'], h['eff_Pt0'],1.,1.,"B")
h['eff_Eta2'].Divide(h['eff_Eta1'],h['eff_Eta0'],1.,1.,"B")
h['eff_P2'].Draw("E")
h['eff_P2'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2'].Draw("E")
h['eff_Pt2'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2'].Draw("E")
h['eff_Eta2'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1_e'].Sumw2()
h['eff_Pt1_e'].Sumw2()
h['eff_Eta1_e'].Sumw2()
h['eff_P2_e'].Divide(h['eff_P1_e'], h['eff_P0_e'],1.,1.,"B")
h['eff_Pt2_e'].Divide(h['eff_Pt1_e'], h['eff_Pt0_e'],1.,1.,"B")
h['eff_Eta2_e'].Divide(h['eff_Eta1_e'],h['eff_Eta0_e'],1.,1.,"B")
h['eff_P2_e'].Draw("E")
h['eff_P2_e'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_e'].Draw("E")
h['eff_Pt2_e'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_e'].Draw("E")
h['eff_Eta2_e'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1_mu'].Sumw2()
h['eff_Pt1_mu'].Sumw2()
h['eff_Eta1_mu'].Sumw2()
h['eff_P2_mu'].Divide(h['eff_P1_mu'], h['eff_P0_mu'],1.,1.,"B")
h['eff_Pt2_mu'].Divide(h['eff_Pt1_mu'], h['eff_Pt0_mu'],1.,1.,"B")
h['eff_Eta2_mu'].Divide(h['eff_Eta1_mu'],h['eff_Eta0_mu'],1.,1.,"B")
h['eff_P2_mu'].Draw("E")
h['eff_P2_mu'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_mu'].Draw("E")
h['eff_Pt2_mu'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_mu'].Draw("E")
h['eff_Eta2_mu'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1_tau'].Sumw2()
h['eff_Pt1_tau'].Sumw2()
h['eff_Eta1_tau'].Sumw2()
h['eff_P2_tau'].Divide(h['eff_P1_tau'], h['eff_P0_tau'],1.,1.,"B")
h['eff_Pt2_tau'].Divide(h['eff_Pt1_tau'], h['eff_Pt0_tau'],1.,1.,"B")
h['eff_Eta2_tau'].Divide(h['eff_Eta1_tau'],h['eff_Eta0_tau'],1.,1.,"B")
h['eff_P2_tau'].Draw("E")
h['eff_P2_tau'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_tau'].Draw("E")
h['eff_Pt2_tau'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_tau'].Draw("E")
h['eff_Eta2_tau'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1_pi'].Sumw2()
h['eff_Pt1_pi'].Sumw2()
h['eff_Eta1_pi'].Sumw2()
h['eff_P2_pi'].Divide(h['eff_P1_pi'], h['eff_P0_pi'],1.,1.,"B")
h['eff_Pt2_pi'].Divide(h['eff_Pt1_pi'], h['eff_Pt0_pi'],1.,1.,"B")
h['eff_Eta2_pi'].Divide(h['eff_Eta1_pi'],h['eff_Eta0_pi'],1.,1.,"B")
h['eff_P2_pi'].Draw("E")
h['eff_P2_pi'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_pi'].Draw("E")
h['eff_Pt2_pi'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_pi'].Draw("E")
h['eff_Eta2_pi'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1_ka'].Sumw2()
h['eff_Pt1_ka'].Sumw2()
h['eff_Eta1_ka'].Sumw2()
h['eff_P2_ka'].Divide(h['eff_P1_ka'], h['eff_P0_ka'],1.,1.,"B")
h['eff_Pt2_ka'].Divide(h['eff_Pt1_ka'], h['eff_Pt0_ka'],1.,1.,"B")
h['eff_Eta2_ka'].Divide(h['eff_Eta1_ka'],h['eff_Eta0_ka'],1.,1.,"B")
h['eff_P2_ka'].Draw("E")
h['eff_P2_ka'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_ka'].Draw("E")
h['eff_Pt2_ka'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_ka'].Draw("E")
h['eff_Eta2_ka'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1_oth'].Sumw2()
h['eff_Pt1_oth'].Sumw2()
h['eff_Eta1_oth'].Sumw2()
h['eff_P2_oth'].Divide(h['eff_P1_oth'], h['eff_P0_oth'],1.,1.,"B")
h['eff_Pt2_oth'].Divide(h['eff_Pt1_oth'], h['eff_Pt0_oth'],1.,1.,"B")
h['eff_Eta2_oth'].Divide(h['eff_Eta1_oth'],h['eff_Eta0_oth'],1.,1.,"B")
h['eff_P2_oth'].Draw("E")
h['eff_P2_oth'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_oth'].Draw("E")
h['eff_Pt2_oth'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_oth'].Draw("E")
h['eff_Eta2_oth'].SetTitle("PseudoRapidity ;#eta;Eff.")

h['eff_P1_mix'].Sumw2()
h['eff_Pt1_mix'].Sumw2()
h['eff_Eta1_mix'].Sumw2()
h['eff_P2_mix'].Divide(h['eff_P1_mix'], h['eff_P0_mix'],1.,1.,"B")
h['eff_Pt2_mix'].Divide(h['eff_Pt1_mix'], h['eff_Pt0_mix'],1.,1.,"B")
h['eff_Eta2_mix'].Divide(h['eff_Eta1_mix'],h['eff_Eta0_mix'],1.,1.,"B")
h['eff_P2_mix'].Draw("E")
h['eff_P2_mix'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_mix'].Draw("E")
h['eff_Pt2_mix'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_mix'].Draw("E")
h['eff_Eta2_mix'].SetTitle("PseudoRapidity ;#eta;Eff.")

r.gStyle.SetOptStat(0)
hfile =tmp1+"_ana.root"
print hfile
r.gROOT.cd()#gr.
ut.writeHists(h,hfile)
