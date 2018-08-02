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

ut.bookHist(h,'DP_pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_pi_oth','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_pi_oth','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_pi_oth','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_ka_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_ka_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_ka_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_ka_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_ka_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_ka_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_2pi0_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_2pi0_oth','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_3pi_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_3pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_4pi_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_4pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_2pi0_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_2pi0_oth','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_3pi_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_3pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_4pi_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_4pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_2pi0_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_2pi0_oth','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_3pi_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_3pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_4pi_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_4pi_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_e_oth','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_e_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_mu_oth','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_mu_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_e_oth','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_e_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_mu_oth','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_mu_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_e_oth','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_e_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_mu_oth','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_mu_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_tau','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_mix','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_mix','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_single','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_single','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_tau','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_mix','invariant Mass(GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_mix','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_single','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_single','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_tau','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_mix','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_mix','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_single','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_single','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_2pi0','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_2pi0','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_3pi','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_3pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_4pi','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_4pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_2pi0','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_2pi0','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_3pi','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_3pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_4pi','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_4pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_2pi0','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_2pi0','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_3pi','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_3pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_4pi','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_4pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)

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
ut.bookHist(h,'eff_P0_e_oth','',50,0.,250.)
ut.bookHist(h,'eff_P1_e_oth','',50,0.,250.)
ut.bookHist(h,'eff_P2_e_oth','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_e_oth','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_e_oth','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_e_oth','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_e_oth','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_e_oth','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_e_oth','',50,2.,8.)
ut.bookHist(h,'eff_P0_mu_oth','',50,0.,250.)
ut.bookHist(h,'eff_P1_mu_oth','',50,0.,250.)
ut.bookHist(h,'eff_P2_mu_oth','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_mu_oth','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_mu_oth','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_mu_oth','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_mu_oth','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_mu_oth','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_mu_oth','',50,2.,8.)
ut.bookHist(h,'eff_P0_pi0','',50,0.,250.)
ut.bookHist(h,'eff_P1_pi0','',50,0.,250.)
ut.bookHist(h,'eff_P2_pi0','',50,0.,250.)
ut.bookHist(h,'eff_Pt0_pi0','',50,0.,8.)
ut.bookHist(h,'eff_Pt1_pi0','',50,0.,8.)
ut.bookHist(h,'eff_Pt2_pi0','',50,0.,8.)
ut.bookHist(h,'eff_Eta0_pi0','',50,2.,8.)
ut.bookHist(h,'eff_Eta1_pi0','',50,2.,8.)
ut.bookHist(h,'eff_Eta2_pi0','',50,2.,8.)


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
    theta=m.atan(r.TMath.Sqrt(m_x**2.+m_y**2.)/m_z)
    return a_x,a_y,theta
def invariant(S_x,S_y,e_0,e_1,p_0,p_1):
    th=r.TMath.Sqrt((S_x[0]-S_x[1])**2.+(S_y[0]-S_y[1])**2.)
    return r.TMath.Sqrt((e_0+e_1)**2.-(p_0**2.+p_1**2.+2.*p_0*p_1*m.cos(th)))


def myEventLoop(n,xsw):
    neg,pos=0,0
    dp,     e,      mu,     tau,    pi,     ka,     oth,    pi0,    oth_d   =0,0,0, 0, 0, 0, 0, 0, 0
    dp_a,   e_a,    mu_a,   tau_a,  pi_a,   ka_a,   oth_a,  pi0_v,  oth_da  =0,0,0, 0, 0, 0, 0, 0, 0
    dp_v,   e_v,    mu_v,   tau_v,  pi_v,   ka_v,   oth_v,  pi0_a,  oth_dv  =0,0,0, 0, 0, 0, 0, 0, 0
    rc=sTree.GetEntry(n)
    had,vec,oth_mo=[],[],[]
    mo=-99
    moth_q=0
    moth_dau=0
    E,P,M,S_y,S_x=[],[],[],[],[]
    if pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,'meson',sTree.MCTrack[0].GetPdgCode())
    wg = sTree.MCTrack[1].GetWeight()
    if abs(wg)<0.00000001: wg = wgnew(sTree.MCTrack[1])
    for mc,track in enumerate(sTree.MCTrack):
        if track.GetMotherId()==1 and (sTree.MCTrack[1].GetPdgCode()==9900015 or sTree.MCTrack[1].GetPdgCode()==4900023):
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
            vtx=r.TVector3(X,Y,Z)
            mom=r.TLorentzVector()
            track.Get4Momentum(mom)
            S_x.append(a[0])
            S_y.append(a[1])
            dp+=1
            if isInFiducial(X,Y,Z):
                dp_v+=1
                h['eff_P0'].Fill(track.GetP())          
                h['eff_Pt0'].Fill(track.GetPt())
                h['eff_Eta0'].Fill(track.GetRapidity())
                if checkFiducialVolume(vtx,a):
                    dp_a+=1
                    h['eff_P1'].Fill(track.GetP())          
                    h['eff_Pt1'].Fill(track.GetPt())
                    h['eff_Eta1'].Fill(track.GetRapidity())
            if abs(track.GetPdgCode())==11 or abs(track.GetPdgCode())==13 or abs(track.GetPdgCode())==15:
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
        if moth_q<4:
        #if len(had)!=0 and moth_q<4:
            for daug in had:
                if track.GetMotherId()==daug:
                    moth_q+=1 
                    X=track.GetStartX()
                    Y=track.GetStartY()
                    Z=track.GetStartZ()
                    Px=track.GetPx()
                    Py=track.GetPy()
                    Pz=track.GetPz()
                    vtx=r.TVector3(X,Y,Z)
                    mom=r.TLorentzVector()
                    track.Get4Momentum(mom)
                    if abs(track.GetPdgCode())==211:
                        pi+=1
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
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
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
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
                    elif track.GetPdgCode()==111:
                        pi0+=1
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        if isInFiducial(X,Y,Z):
                             pi0_v+=1
                             h['eff_P0_pi0'].Fill(track.GetP(),wg)            
                             h['eff_Pt0_pi0'].Fill(track.GetPt(),wg)
                             h['eff_Eta0_pi0'].Fill(track.GetRapidity(),wg)
                             if ang_hcal(vtx,a):
                                 pi0_a+=1
                                 h['eff_P1_pi0'].Fill(track.GetP())          
                                 h['eff_Pt1_pi0'].Fill(track.GetPt())
                                 h['eff_Eta1_pi0'].Fill(track.GetRapidity())

                    else:
                        oth_mo.append(mc)
                        oth+=1
                        vec.append(track.GetPdgCode())
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        if isInFiducial(X,Y,Z):
                            oth_v+=1
                            h['eff_P0_oth'].Fill(track.GetP())          
                            h['eff_Pt0_oth'].Fill(track.GetPt())
                            h['eff_Eta0_oth'].Fill(track.GetRapidity())
                            if checkFiducialVolume(vtx,a):
                                oth_a+=1
                                h['eff_P1_oth'].Fill(track.GetP())          
                                h['eff_Pt1_oth'].Fill(track.GetPt())
                                h['eff_Eta1_oth'].Fill(track.GetRapidity())
        if moth_dau<8:
        #if len(oth_mo)!=0 and moth_dau<8:
            for daug in oth_mo:
                if track.GetMotherId()==daug:
                    moth_dau+=1
                    X=track.GetStartX()
                    Y=track.GetStartY()
                    Z=track.GetStartZ()
                    Px=track.GetPx()
                    Py=track.GetPy()
                    Pz=track.GetPz()
                    vtx=r.TVector3(X,Y,Z)
                    mom=r.TLorentzVector()
                    track.Get4Momentum(mom)
                    if abs(track.GetPdgCode())==211:
                        pi+=1
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        if isInFiducial(X,Y,Z):
                            pi_v+=1
                            if ang_ecal(vtx,a):
                                pi_a+=1
                    elif abs(track.GetPdgCode())==321:
                        ka+=1
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        if isInFiducial(X,Y,Z):
                            ka_v+=1
                            if ang_ecal(vtx,a):
                                ka_a+=1
                    elif track.GetPdgCode()==111:
                        pi0+=1
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        if isInFiducial(X,Y,Z):
                             pi0_v+=1
                             if ang_hcal(vtx,a):
                                 pi0_a+=1
                    elif abs(track.GetPdgCode())==11:
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        e+=1
                        if isInFiducial(X,Y,Z):
                            e_v+=1
                            if ang_ecal(vtx,a):
                                e_a+=1
                    elif abs(track.GetPdgCode())==13:
                        mu+=1
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        if isInFiducial(X,Y,Z):
                            mu_v+=1
                            if ang_muon(vtx,a):
                                mu_a+=1
                    else:
                        oth_d+=1
                        vec.append(track.GetPdgCode())
                        P.append(track.GetP())
                        E.append(track.GetEnergy())
                        M.append(mom.M())
                        a=angle(track.GetP(),Px,Py,Pz)
                        S_x.append(a[0])
                        S_y.append(a[1])
                        if isInFiducial(X,Y,Z):
                            oth_dv+=1
                            if checkFiducialVolume(vtx,a):
                                oth_da+=1

    
    if pi0==0 and ka==0 and (oth==0 or oth_d==0) and e==2 and mu==0 and tau==0 and pi==0:  
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_e'].Fill(mass)
        h['DPW_e'].Fill(mass,wg*xsw)
        if e_v==2:  
            h['DPves_e'].Fill(mass)
            h['DPvesW_e'].Fill(mass,wg*xsw)
            if e_a==2:
                h['DPang_e'].Fill(mass)
                h['DPangW_e'].Fill(mass,wg*xsw)

    elif pi0==0 and ka==0 and (oth==0 or oth_d==0) and e==0 and mu==2 and tau==0 and pi==0:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_mu'].Fill(mass)
        h['DPW_mu'].Fill(mass,wg*xsw)
        if mu_v==2:
            h['DPves_mu'].Fill(mass)
            h['DPvesW_mu'].Fill(mass,wg*xsw)
            if mu_a==2:
                h['DPang_mu'].Fill(mass)
                h['DPangW_mu'].Fill(mass,wg*xsw)

    elif pi0==0 and ka==0 and (oth==0 or oth_d==0) and e==0 and mu==0 and tau==2 and pi==0:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_tau'].Fill(mass)
        h['DPW_tau'].Fill(mass,wg*xsw)
        if tau_v==2:
            h['DPves_tau'].Fill(mass)
            h['DPvesW_tau'].Fill(mass,wg*xsw)
            if tau_a==2:
                h['DPang_tau'].Fill(mass)
                h['DPangW_tau'].Fill(mass,wg*xsw)

    elif pi0==0 and ka==0 and (oth==0 or oth_d==0) and e==0 and mu==0 and tau==0 and pi==2:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_pi'].Fill(mass)
        h['DPW_pi'].Fill(mass,wg*xsw)
        if pi_v==2:
            h['DPves_pi'].Fill(mass)
            h['DPvesW_pi'].Fill(mass,wg*xsw)
            if pi_a==2:
                h['DPang_pi'].Fill(mass)
                h['DPangW_pi'].Fill(mass,wg*xsw)

    elif (oth==0 or oth_d==0) and pi0==0 and pi==0 and e==0 and mu==0 and tau==0 and ka==2:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_ka'].Fill(mass)       
        h['DPW_ka'].Fill(mass,wg*xsw)
        if ka_v==2:
            h['DPves_ka'].Fill(mass)        
            h['DPvesW_ka'].Fill(mass,wg*xsw)
            if ka_a==2:
                h['DPang_ka'].Fill(mass)        
                h['DPangW_ka'].Fill(mass,wg*xsw)

    elif oth==2 and oth_d!=0 and pi==0 and pi0==0 and e==0 and mu==0 and tau==0 and ka==0:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        if abs(vec[0])==abs(vec[1]):
            h['DP_oth'].Fill(mass)       
            h['DPW_oth'].Fill(mass,wg*xsw)
        else:
            h['DP_mix'].Fill(mass)       
            h['DPW_mix'].Fill(mass,wg*xsw)
        if oth_v==2: 
            if abs(vec[0])==abs(vec[1]):
                h['DPves_oth'].Fill(mass)       
                h['DPvesW_oth'].Fill(mass,wg*xsw)
            else:
                h['DPves_mix'].Fill(mass)       
                h['DPvesW_mix'].Fill(mass,wg*xsw)
            if oth_a==2:
                if abs(vec[0])==abs(vec[1]):
                    h['DPang_oth'].Fill(mass)       
                    h['DPangW_oth'].Fill(mass,wg*xsw)
                else:
                    h['DPang_mix'].Fill(mass)       
                    h['DPangW_mix'].Fill(mass,wg*xsw)

    elif pi==2 and pi0==2 and (oth==0 or oth_d==0) and ka==0 and e==0 and mu==0 and tau==0:
        #mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_2pi0'].Fill(mass_mc)       
        h['DPW_2pi0'].Fill(mass_mc,wg*xsw)
        if pi_v==2 and pi0_v==2:
            h['DPves_2pi0'].Fill(mass_mc)        
            h['DPvesW_2pi0'].Fill(mass_mc,wg*xsw)
            if pi_a==2 and pi0_a==2:
                h['DPang_2pi0'].Fill(mass_mc)        
                h['DPangW_2pi0'].Fill(mass_mc,wg*xsw)

    elif pi==2 and pi0==1 and (oth==0 or oth_d==0) and ka==0 and e==0 and mu==0 and tau==0:
        #mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_3pi'].Fill(mass_mc)       
        h['DPW_3pi'].Fill(mass_mc,wg*xsw)
        if pi_v==2 and pi0_v==1:
            h['DPves_3pi'].Fill(mass_mc)        
            h['DPvesW_3pi'].Fill(mass_mc,wg*xsw)
            if pi_a==2 and pi0_a==1:
                h['DPang_3pi'].Fill(mass_mc)        
                h['DPangW_3pi'].Fill(mass_mc,wg*xsw)

    elif pi==4 and pi0==0 and e==0 and mu==0 and tau==0 and (oth==0 or oth_d==0) and ka==0:
        #mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['DP_4pi'].Fill(mass_mc)       
        h['DPW_4pi'].Fill(mass_mc,wg*xsw)
        if pi_v==4:
            h['DPves_4pi'].Fill(mass_mc)        
            h['DPvesW_4pi'].Fill(mass_mc,wg*xsw)
            if pi_a==4:
                h['DPang_4pi'].Fill(mass_mc)        
                h['DPangW_4pi'].Fill(mass_mc,wg*xsw)
 
    else:
        #print n,pi,ka,oth
        h['DP_single'].Fill(mass_mc)        
        h['DPW_single'].Fill(mass_mc,wg*xsw)
        if oth_v!=2 and oth_v!=0:
            h['DPves_single'].Fill(mass_mc)     
            h['DPvesW_single'].Fill(mass_mc,wg*xsw)    
            if oth_a!=2 and oth_a!=0:
                h['DPang_single'].Fill(mass_mc)     
                h['DPangW_single'].Fill(mass_mc,wg*xsw)

    if dp==2:
        mass=invariant(S_x, S_y,E[0],E[1],P[0],P[1])
        h['eff_Pw'].Fill(sTree.MCTrack[1].GetP(),wg)            
        h['eff_Ptw'].Fill(sTree.MCTrack[1].GetPt(),wg)
        h['eff_Etaw'].Fill(sTree.MCTrack[1].GetRapidity(),wg)
        h['DP'].Fill(mass)
        h['DPW'].Fill(mass,wg*xsw)
        if dp_v==2 and ((e_v==2 and mu_v==0 and tau_v==0 and pi_v==0 and ka_v==0 and pi0_v==0 and (oth_v==0 or oth_dv==0)) or (e_v==0 and mu_v==2 and tau_v==0 and pi_v==0 and ka_v==0 and pi0_v==0 and (oth_v==0 or oth_dv==0)) or (e_v==0 and mu_v==0 and tau_v==2 and pi_v==0 and ka_v==0 and pi0_v==0 and (oth_v==0 or oth_dv==0)) or (e_v==0 and mu_v==0 and tau_v==0 and pi_v==2 and ka_v==0 and pi0_v==0 and (oth_v==0 or oth_dv==0)) or (e_v==0 and mu_v==0 and tau_v==0 and pi_v==0 and ka_v==2 and pi0_v==0 and (oth_v==0 or oth_dv==0)) or (e_v==0 and mu_v==0 and tau_v==0 and pi_v==2 and ka_v==0 and pi0_v==1 and (oth_v==0 or oth_dv==0)) or (e_v==0 and mu_v==0 and tau_v==0 and pi_v==2 and ka_v==0 and pi0_v==2 and (oth_v==0 or oth_dv==0)) or (e_v==0 and mu_v==0 and tau_v==0 and pi_v==4 and ka_v==0 and pi0_v==0 and (oth_v==0 or oth_dv==0))):
            h['DPves'].Fill(mass)
            h['DPvesW'].Fill(mass,wg*xsw)
            h['eff_P0w'].Fill(sTree.MCTrack[1].GetP(),wg*xsw)           
            h['eff_Pt0w'].Fill(sTree.MCTrack[1].GetPt(),wg*xsw)
            h['eff_Eta0w'].Fill(sTree.MCTrack[1].GetRapidity(),wg*xsw)
            if dp_a==2 and ((e_a==2 and mu_a==0 and tau_a==0 and pi_a==0 and ka_a==0 and pi0_a==0 and (oth_a==0 or oth_da==0)) or (e_a==0 and mu_a==2 and tau_a==0 and pi_a==0 and ka_a==0 and pi0_a==0 and (oth_a==0 or oth_da==0)) or (e_a==0 and mu_a==0 and tau_a==2 and pi_a==0 and ka_a==0 and pi0_a==0 and (oth_a==0 or oth_da==0)) or (e_a==0 and mu_a==0 and tau_a==0 and pi_a==2 and ka_a==0 and pi0_a==0 and (oth_a==0 or oth_da==0)) or (e_a==0 and mu_a==0 and tau_a==0 and pi_a==0 and ka_a==2 and pi0_a==0 and (oth_a==0 or oth_da==0)) or (e_a==0 and mu_a==0 and tau_a==0 and pi_a==2 and ka_a==0 and pi0_a==1 and (oth_a==0 or oth_da==0)) or (e_a==0 and mu_a==0 and tau_a==0 and pi_a==2 and ka_a==0 and pi0_a==2 and (oth_a==0 or oth_da==0)) or (e_a==0 and mu_a==0 and tau_a==0 and pi_a==4 and ka_a==0 and pi0_a==0 and (oth_a==0 or oth_da==0))):
                h['DPang'].Fill(mass)
                h['DPangW'].Fill(mass,wg*xsw)
                h['eff_P1w'].Fill(sTree.MCTrack[1].GetP(),wg*xsw)           
                h['eff_Pt1w'].Fill(sTree.MCTrack[1].GetPt(),wg*xsw)
                h['eff_Eta1w'].Fill(sTree.MCTrack[1].GetRapidity(),wg*xsw)

nEvents =sTree.GetEntries()
print nEvents
for n in range(nEvents):
    if not pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,pro,0)
    myEventLoop(n,xsw)

o1 = tmp1+"_true.dat"
o2 = tmp1+"_e.dat"
o3 = tmp1+"_mu.dat"
o4 = tmp1+"_pi.dat"
o5 = tmp1+"_ka.dat"
o6 = tmp1+"_oth.dat"
o7 = tmp1+"_mix.dat"
o8 = tmp1+"_tau.dat"
#o9 = tmp1+"_single.dat"
o24 = tmp1+"_4pi.dat"
o25 = tmp1+"_3pi.dat"
o26 = tmp1+"_2pi0.dat"
o27 = tmp1+"_e_oth.dat"
o28 = tmp1+"_mu_oth.dat"
a=open(o1,'w+')
b=open(o2,'w+')
c=open(o3,'w+')
d=open(o4,'w+')
e=open(o5,'w+')
f=open(o6,'w+')
g=open(o7,'w+')
t=open(o8,'w+')
#s=open(o9,w+)
d2=open(o24,'w+')
e2=open(o25,'w+')
f2=open(o26,'w+')
o11 = tmp1+"_true_table.dat"
o12 = tmp1+"_e_table.dat"
o13 = tmp1+"_mu_table.dat"
o14 = tmp1+"_pi_table.dat"
o15 = tmp1+"_ka_table.dat"
o16 = tmp1+"_oth_table.dat"
o17 = tmp1+"_mix_table.dat"
o18 = tmp1+"_tau_table.dat"
o19 = tmp1+"_single_table.dat"
o124 = tmp1+"_4pi_table.dat"
o125 = tmp1+"_3pi_table.dat"
o126 = tmp1+"_2pi0_table.dat"
a1=open(o11,'w+')
b1=open(o12,'w+')
c1=open(o13,'w+')
d1=open(o14,'w+')
e1=open(o15,'w+')
f1=open(o16,'w+')
g1=open(o17,'w+')
t1=open(o18,'w+')
s1=open(o19,'w+')
d12=open(o124,'w+')
e12=open(o125,'w+')
f12=open(o126,'w+')

if float(h['DPvesW'].Integral())!=0: Acc=float(h['DPangW'].Integral())/float(h['DPvesW'].Integral())
if float(h['DPvesW_e'].Integral())!=0: Acc_e=float(h['DPangW_e'].Integral())/float(h['DPvesW_e'].Integral())
if float(h['DPvesW_mu'].Integral())!=0: Acc_mu=float(h['DPangW_mu'].Integral())/float(h['DPvesW_mu'].Integral())
if float(h['DPvesW_tau'].Integral())!=0: Acc_tau=float(h['DPangW_tau'].Integral())/float(h['DPvesW_tau'].Integral())
if float(h['DPvesW_pi'].Integral())!=0: Acc_pi=float(h['DPangW_pi'].Integral())/float(h['DPvesW_pi'].Integral())
if float(h['DPvesW_ka'].Integral())!=0: Acc_ka=float(h['DPangW_ka'].Integral())/float(h['DPvesW_ka'].Integral())
if float(h['DPvesW_oth'].Integral())!=0: Acc_oth=float(h['DPangW_oth'].Integral())/float(h['DPvesW_oth'].Integral())
if float(h['DPvesW_mix'].Integral())!=0: Acc_mix=float(h['DPangW_mix'].Integral())/float(h['DPvesW_mix'].Integral())
if float(h['DPvesW_2pi0'].Integral())!=0: Acc_2pi0=float(h['DPangW_2pi0'].Integral())/float(h['DPvesW_2pi0'].Integral())
if float(h['DPvesW_3pi'].Integral())!=0: Acc_3pi=float(h['DPangW_3pi'].Integral())/float(h['DPvesW_3pi'].Integral())
if float(h['DPvesW_4pi'].Integral())!=0: Acc_4pi=float(h['DPangW_4pi'].Integral())/float(h['DPvesW_4pi'].Integral())
if float(h['DPvesW_e_oth'].Integral())!=0: Acc_e_oth=float(h['DPangW_e_oth'].Integral())/float(h['DPvesW_e_oth'].Integral())
if float(h['DPvesW_mu_oth'].Integral())!=0: Acc_mu_oth=float(h['DPangW_mu_oth'].Integral())/float(h['DPvesW_mu_oth'].Integral())


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
if float(h['DPves_2pi0'].Integral())!=0:
    RecW_2pi0=float(h['DPangW_2pi0'].Integral())/float(h['DPves_2pi0'].Integral())*2e20
    f2.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_2pi0, RecW_2pi0, Acc_2pi0))
    f2.write('\n')
if float(h['DPves_3pi'].Integral())!=0:
    RecW_3pi=float(h['DPangW_3pi'].Integral())/float(h['DPves_3pi'].Integral())*2e20
    e2.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_3pi, RecW_3pi, Acc_3pi)) 
    e2.write('\n')
if float(h['DPves_4pi'].Integral())!=0:
    RecW_4pi=float(h['DPangW_4pi'].Integral())/float(h['DPves_4pi'].Integral())*2e20
    d2.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_4pi, RecW_4pi, Acc_4pi))
    d2.write('\n')


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
s1.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_single'].Integral()), float(h['DPves_single'].Integral()), float(h['DPang_single'].Integral())))
s1.write('\n')
f12.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_2pi0'].Integral()), float(h['DPves_2pi0'].Integral()), float(h['DPang_2pi0'].Integral())))
f12.write('\n')
e12.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_3pi'].Integral()), float(h['DPves_3pi'].Integral()), float(h['DPang_3pi'].Integral())))
e12.write('\n')
d12.write('%.4g %s %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP_4pi'].Integral()), float(h['DPves_4pi'].Integral()), float(h['DPang_4pi'].Integral())))
d12.write('\n')

a.close()
b.close()
c.close()
d.close()
e.close()
f.close()
t.close()
g.close()
#s.close()
d2.close()
e2.close()
f2.close()
a1.close()
b1.close()
c1.close()
d1.close()
e1.close()
f1.close()
t1.close()
g1.close()
s1.close()
d12.close()
e12.close()
f12.close()

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
h['eff_P1_pi0'].Sumw2()
h['eff_Pt1_pi0'].Sumw2()
h['eff_Eta1_pi0'].Sumw2()
h['eff_P2_pi0'].Divide(h['eff_P1_pi0'], h['eff_P0_pi0'],1.,1.,"B")
h['eff_Pt2_pi0'].Divide(h['eff_Pt1_pi0'], h['eff_Pt0_pi0'],1.,1.,"B")
h['eff_Eta2_pi0'].Divide(h['eff_Eta1_pi0'],h['eff_Eta0_pi0'],1.,1.,"B")
h['eff_P2_pi0'].Draw("E")
h['eff_P2_pi0'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2_pi0'].Draw("E")
h['eff_Pt2_pi0'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2_pi0'].Draw("E")
h['eff_Eta2_pi0'].SetTitle("PseudoRapidity ;#eta;Eff.")


r.gStyle.SetOptStat(0)
hfile =tmp1+"_ana.root"
print hfile
r.gROOT.cd()#gr.
ut.writeHists(h,hfile)
