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
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:g:", ["nEvents=","geoFile="])
except getopt.GetoptError:
    print 'no file'
    sys.exit()
for o,a in opts:
    if o in ('-f',): inputFile = a
    if o in ('-g', '--geoFile',): geoFile = a
 
tmp=inputFile.replace("/eos/experiment/ship/data/DarkPhoton/PBC-June-3/rec/","")
tmp1=tmp.replace('_rec.root','')
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
ecalGeoFile = ShipGeo.ecal.File
hcalGeoFile = ShipGeo.hcal.File 
dy = ShipGeo.Yheight/u.m
MeasCut=25


# -----Create geometry----------------------------------------------
import shipDet_conf
run = r.FairRunSim()
run.SetName("TGeant4")  # Transport engine
run.SetOutputFile(ROOT.TMemFile('output', 'recreate'))  # Output file
run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for the mag field
rtdb = run.GetRuntimeDb()
# -----Create geometry----------------------------------------------
modules = shipDet_conf.configure(run,ShipGeo)

import geomGeant4
if hasattr(ShipGeo.Bfield,"fieldMap"):
  fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True, withVirtualMC = False)
else:
  print "no fieldmap given, geofile too old, not anymore support"
  exit(-1)
sGeo   = fgeo.FAIRGeom
geoMat =  r.genfit.TGeoMaterialInterface()
ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
bfield = r.genfit.FairShipFields()
bfield.setField(fieldMaker.getGlobalField())
fM = r.genfit.FieldManager.getInstance()
fM.init(bfield)

volDict = {}
i=0
for x in ROOT.gGeoManager.GetListOfVolumes():
 volDict[i]=x.GetName()
 i+=1

# prepare veto decisions
import shipVeto
veto = shipVeto.Task(sTree)
vetoDets={}

PDG = r.TDatabasePDG.Instance()
 
h={}
ut.bookHist(h,'DOCA','DOCA')
ut.bookHist(h,'DP','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DP_e','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW_e','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves_e','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW_e','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang_e','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW_e','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DP_mu','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW_mu','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves_mu','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW_mu','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang_mu','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW_mu','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DP_tau','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW_tau','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves_tau','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW_tau','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang_tau','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW_tau','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DP_pi','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW_pi','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves_pi','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW_pi','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang_pi','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW_pi','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DP_ka','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW_ka','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves_ka','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW_ka','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang_ka','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW_ka','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DP_oth','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW_oth','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves_oth','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW_oth','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang_oth','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW_oth','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DP_mix','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPW_mix','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPves_mix','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPvesW_mix','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPang_mix','invariant Mass (GeV)',50,0.,mass_mc+5.)
ut.bookHist(h,'DPangW_mix','invariant Mass with Weights (GeV)',50,0.,mass_mc+5.)

 
import TrackExtrapolateTool 
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
def checkFiducialVolume(sTree,tkey,dy):
# extrapolate track to middle of magnet and check if in decay volume
   inside = True
   #if not fiducialCut: return True
   fT = sTree.FitTracks[tkey]
   rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,ShipGeo.Bfield.z)
   if not rc: return False
   if not dist2InnerWall(pos.X(),pos.Y(),pos.Z())>0: return False
   return inside

def isInFiducial(X,Y,Z):
    if Z > ShipGeo.TrackStation1.z : return False
    if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
    if dist2InnerWall(X,Y,Z)<5*u.cm: return False
    return True 
 
def findmum():
    for dp_ind,dp_tr in enumerate(sTree.MCTrack):
        if dp_tr.GetPdgCode()==9900015 or dp_tr.GetPdgCode()==4900023:
            mum_id=dp_tr.GetMotherId()
            dp_id=dp_ind
            #print dp_id
            if pro=='qcd' and dp_id==0: continue
            #print mum_id 
            mum_pdg=sTree.MCTrack[mum_id].GetPdgCode()
            if pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,'meson',mum_pdg)
            else: xsw = dputil.getDPprodRate(mass_mc,eps,pro,0) 
            wg = sTree.MCTrack[dp_id].GetWeight()
            #print dp_id 
            dp_mom=r.TVector3(sTree.MCTrack[dp_id].GetPx(),sTree.MCTrack[dp_id].GetPy(),sTree.MCTrack[dp_id].GetPz())
            dp_mag=sTree.MCTrack[dp_id].GetP()
            break
        else: xsw,wg,dp_id,dp_mom,dp_mag=-99,-99,-99,-99,-99
    return xsw,wg,dp_id,dp_mom,dp_mag

def ImpactParameter(point,tPos,tMom):
    t = 0
    if hasattr(tMom,'P'): P = tMom.P()
    else:                 P = tMom.Mag()
    for i in range(3):   t += tMom(i)/P*(point(i)-tPos(i))
    dist = 0
    for i in range(3):   dist += (point(i)-tPos(i)-t*tMom(i)/P)**2
    dist = r.TMath.Sqrt(dist)
    return dist

def IP(X,Y,Z,dp_M,dp_Mag):
    target = r.TVector3(0., 0., ShipGeo.target.z0)
    delta = 0.
    delta += (target(0) - X) * dp_M(0)/dp_Mag
    delta += (target(1) - Y) * dp_M(1)/dp_Mag
    delta += (target(2) - Z) * dp_M(2)/dp_Mag
    ip = 0.
    ip += (target(0) - X - delta*dp_M(0)/dp_Mag)**2.
    ip += (target(1) - Y - delta*dp_M(1)/dp_Mag)**2.
    ip += (target(2) - Z - delta*dp_M(2)/dp_Mag)**2.
    return r.TMath.Sqrt(ip)

def checkOrigin(tr_MOM, dp_id, mc):
    Lo=0
    while Lo<11:
        pid=sTree.MCTrack[mc].GetPdgCode()
        mc=tr_MOM
        if tr_MOM == int(dp_id):
            return True, pid, sTree.MCTrack[tr_MOM].GetPdgCode()
        elif abs(sTree.MCTrack[tr_MOM].GetPdgCode())<7 or abs(sTree.MCTrack[tr_MOM].GetPdgCode())==21:
            return True, pid, sTree.MCTrack[tr_MOM].GetPdgCode()
        else: 
            Lo+=1
            tr_MOM=sTree.MCTrack[tr_MOM].GetMotherId()
    if Lo>10: return False, -99, -99
def checkTrue(sTree, dp_id):
    PID=[]
    Lo, Ko = 0,0
    for mc,tr in enumerate(sTree.MCTrack):
        pid = tr.GetPdgCode()
        mom = tr.GetMotherId()
        if mc>1:
            mom_pid=sTree.MCTrack[mom].GetPdgCode()
            if abs(pid)>9 and abs(pid)!=21:
                if mom==dp_id:#leptons
                    Lo+=1
                    PID.append(pid)
                    #print n, "th event; Daughters ",pid
                elif (abs(pid)==211 or abs(pid)==111 or abs(pid)==321 ) and (abs(mom_pid)<7 or abs(mom_pid)==21):#hadrons
                    Ko+=1
                    PID.append(pid)
                    #print n, "th event; Daughters ",pid
        #if Lo==2: return PID
    return PID
def myEventLoop(n):
    rc=sTree.GetEntry(n)
    fm=findmum()
    xsw=fm[0]
    wg=fm[1]
    dp_id=fm[2]
    dp_M=fm[3]
    dp_Mag=fm[4]
    MA,MAS=[],[] 
    DPmom=r.TLorentzVector(0.,0.,0.,0.)
    DPma=r.TLorentzVector(0.,0.,0.,0.) 
    dau=-99
    DOC=[]
    TR=0
    T1,T2=[],[]
    VES=0 
    r_track=0
    f_track=0
    RECO=0
    debug=0
    e_a, mu_a, tau_a, pi_a, pi0_a, ka_a, oth_a = 0, 0, 0, 0, 0, 0, 0
    e_v, mu_v, tau_v, pi_v, pi0_v, ka_v, oth_v = 0, 0, 0, 0, 0, 0, 0
    e,   mu,   tau,   pi,   pi0,   ka,   oth   = 0, 0, 0, 0, 0, 0, 0
    if xsw==-99 and wg==-99 and dp_id==-99: return
    h['DP'].Fill(mass_mc) 
    dau=checkTrue(sTree,dp_id)
    for xxx in dau:
        TR+=1
        if abs(xxx)==11:
            e+=1
        elif abs(xxx)==13:
            mu+=1
        elif abs(xxx)==15:
            tau+=1
        elif abs(xxx)==211:
            pi+=1
        elif abs(xxx)==111:
            pi0+=1
        elif abs(xxx)==321:
            ka+=1
        else:
            oth+=1
            print n, xxx
    for F,FIT in enumerate(sTree.FitTracks):
        fitStatus = FIT.getFitStatus()
        if not fitStatus.isFitConverged(): continue
        xx = FIT.getFittedState()
        mc = sTree.MCTrack[sTree.fitTrack2MC[F]]
        vtx=r.TVector3(mc.GetStartX(), mc.GetStartY(), mc.GetStartZ())
        mom=r.TVector3(mc.GetPx(), mc.GetPy(), mc.GetPz())
        if not isInFiducial(vtx.X(),vtx.Y(),vtx.Z()): continue
        rc, pid, mom_pid= checkOrigin(mc.GetMotherId(),dp_id,sTree.fitTrack2MC[F])
        #f_track+=1
        if abs(sTree.MCTrack[2].GetPdgCode())==15:
            f_track+=1
            tau_v+=1
        if abs(sTree.MCTrack[2].GetPdgCode())!=15:
            if not rc: continue
            if abs(pid)==11:
                f_track+=1
                e_v+=1
            elif abs(pid)==13:
                f_track+=1
                mu_v+=1
            elif abs(pid)==211:
                f_track+=1
                pi_v+=1
            elif abs(pid)==321:
                f_track+=1
                ka_v+=1
            elif abs(pid)==22 or abs(pid)==111:
                f_track+=1
                pi0_v+=1
            else:
                #print xx.getPDG(), sTree.MCTrack[mc.GetMotherId()].GetPdgCode()  
                if abs(xx.getPDG())==11:
                    f_track+=1
                    e_v+=1
                elif abs(xx.getPDG())==13:
                    f_track+=1
                    mu_v+=1
                elif abs(xx.getPDG())==211:
                    f_track+=1
                    pi_v+=1
                elif abs(xx.getPDG())==321:
                    f_track+=1
                    ka_v+=1
                elif abs(xx.getPDG())==22 or abs(xx.getPDG())==111:
                    f_track+=1
                    pi0_v+=1
                #if not checkOrigin(mc.GetMotherId(),dp_id,sTree.fitTrack2MC[F]): continue
                else: oth_v+=1
        if not checkFiducialVolume(sTree,F,dy): continue
        nmeas = fitStatus.getNdf()
        chi2 = fitStatus.getChi2()
        if not nmeas>25: continue
        if not chi2/nmeas<5: continue
        if not xx.getMomMag()>1.: continue
        #RECO+=1
        if abs(sTree.MCTrack[2].GetPdgCode())==15:
            RECO+=1
            tau_a+=1
        if abs(sTree.MCTrack[2].GetPdgCode())!=15:
            #rc, pid, mom_pid= checkOrigin(mc.GetMotherId(),dp_id,sTree.fitTrack2MC[F])
            if not rc: continue
            if abs(pid)==11:
                RECO+=1
                e_a+=1
            elif abs(pid)==13:
                RECO+=1
                mu_a+=1
            elif abs(pid)==211:
                RECO+=1
                pi_a+=1
            elif abs(pid)==321:
                RECO+=1
                ka_a+=1
            elif abs(pid)==22 or abs(pid)==111:
                RECO+=1
                pi0_a+=1
            else:
                #print xx.getPDG(), sTree.MCTrack[mc.GetMotherId()].GetPdgCode()  
                if abs(xx.getPDG())==11:
                    RECO+=1
                    e_a+=1
                elif abs(xx.getPDG())==13:
                    RECO+=1
                    mu_a+=1
                elif abs(xx.getPDG())==211:
                    RECO+=1
                    pi_a+=1
                elif abs(xx.getPDG())==321:
                    RECO+=1
                    ka_a+=1
                elif abs(xx.getPDG())==22 or abs(xx.getPDG())==111:
                    RECO+=1
                    pi0_a+=1
                #if not checkOrigin(mc.GetMotherId(),dp_id,sTree.fitTrack2MC[F]): continue
                else: oth_a+=1
    if TR>1:
        if e>1: h['DP_e'].Fill(mass_mc)
        elif mu>1: h['DP_mu'].Fill(mass_mc)
        elif tau>1: h['DP_tau'].Fill(mass_mc)
        elif pi>1 or (pi0>0 and pi>1): h['DP_pi'].Fill(mass_mc)
        elif ka>1: h['DP_ka'].Fill(mass_mc)
        elif oth>1: h['DP_oth'].Fill(mass_mc)
        else: 
            h['DP_mix'].Fill(mass_mc)
            #print "no event"

    if TR>0:
        #h['DP'].Fill(mass_mc) 
        #print e, mu, tau, pi, pi0, ka, oth 
        if f_track>0:
            h['DPvesW'].Fill(mass_mc,wg)
            h['DPves'].Fill(mass_mc)
            if RECO>0: 
                h['DPangW'].Fill(mass_mc)       
                h['DPang'].Fill(mass_mc,wg*xsw) 
    if f_track>1: 
        if (pi0_v>0 and pi_v>1) or pi_v>1:
            h['DPvesW_pi'].Fill(mass_mc)               
            h['DPves_pi'].Fill(mass_mc,wg) 
            if (pi_a >1 and pi0_a>0) or (pi_a>1):
                h['DPang_pi'].Fill(mass_mc,wg*xsw)
                h['DPangW_pi'].Fill(mass_mc,wg)
            elif (e_a+mu_a+tau_a+pi_a+pi0_a+ka_a)>1: 
                h['DPang_mix'].Fill(mass_mc,wg*xsw) 
                h['DPangW_mix'].Fill(mass_mc,wg)
        elif ka_v>1:
            h['DPves_ka'].Fill(mass_mc,wg)
            h['DPvesW_ka'].Fill(mass_mc)         
            if ka_a>1:
                h['DPang_ka'].Fill(mass_mc,wg*xsw) 
                h['DPangW_ka'].Fill(mass_mc,wg) 
            elif (e_a+mu_a+tau_a+pi_a+pi0_a+ka_a)>1: 
                h['DPang_mix'].Fill(mass_mc,wg*xsw) 
                h['DPangW_mix'].Fill(mass_mc,wg)
        elif tau_v>1:
            h['DPves_tau'].Fill(mass_mc,wg)
            h['DPvesW_tau'].Fill(mass_mc)
            if tau_a>1:
                h['DPang_tau'].Fill(mass_mc,wg*xsw)                
                h['DPangW_tau'].Fill(mass_mc,wg)                     
            elif (e_a+mu_a+tau_a+pi_a+pi0_a+ka_a)>1:
                h['DPang_mix'].Fill(mass_mc,wg*xsw) 
                h['DPangW_mix'].Fill(mass_mc,wg) 
        elif mu_v>1:
            h['DPves_mu'].Fill(mass_mc,wg)
            h['DPvesW_mu'].Fill(mass_mc) 
            if mu_a>1:
                h['DPang_mu'].Fill(mass_mc,wg*xsw)
                h['DPangW_mu'].Fill(mass_mc,wg)
            elif (e_a+mu_a+tau_a+pi_a+pi0_a+ka_a)>1:
                h['DPang_mix'].Fill(mass_mc,wg*xsw) 
                h['DPangW_mix'].Fill(mass_mc,wg)
        elif e_v>1:
            h['DPves_e'].Fill(mass_mc,wg)
            h['DPvesW_e'].Fill(mass_mc)  
            if e_a>1:
                h['DPang_e'].Fill(mass_mc,wg*xsw)
                h['DPangW_e'].Fill(mass_mc,wg)
            elif (e_a+mu_a+tau_a+pi_a+pi0_a+ka_a)>1:
                h['DPang_mix'].Fill(mass_mc,wg*xsw) 
                h['DPangW_mix'].Fill(mass_mc,wg)
        elif (e_v+mu_v+tau_v+pi_v+pi0_v+ka_v)>1:                   
            h['DPvesW_mix'].Fill(mass_mc)                                 
            h['DPves_mix'].Fill(mass_mc,wg)
            if (e_a+mu_a+tau_a+pi_a+pi0_a+ka_a)>1:
                h['DPang_mix'].Fill(mass_mc,wg*xsw) 
                h['DPangW_mix'].Fill(mass_mc,wg)
            #else:
                #print "else_a; ", e_a, mu_a, tau_a, pi_a, pi0_a, ka_a, oth_a
        else:                           
            #print "else_v;  ", e_v, mu_v, tau_v, pi_v, pi0_v, ka_v, oth_v
            h['DPvesW_oth'].Fill(mass_mc)
            #if oth_a<2 and ka_a<2 and (pi_a<2 or pi0_a<2) and mu_a<2 and e_a<2 and tau_a<2:
                #print "else_a; ", e_a, mu_a, tau_a, pi_a, pi0_a, ka_a, oth_a
    
                           
    

nEvents =sTree.GetEntries()
print nEvents
#nEvents = 200
#Dump(sTree.MCTrack)
for n in range(nEvents):
    myEventLoop(n)
 
#z1 = tmp1+"_Ana.dat"
#z=open(z1,'w+')

print "tau_v", h['DPvesW_tau'].Integral()
print "e_v",   h['DPvesW_e'].Integral()
print "mu_v",  h['DPvesW_mu'].Integral()
print "pi_v",  h['DPvesW_pi'].Integral()
print "ka_v",  h['DPvesW_ka'].Integral()
print "mix_v", h['DPvesW_mix'].Integral()
print "oth_v", h['DPvesW_oth'].Integral()
print "tau_a", h['DPangW_tau'].Integral()
print "e_a",   h['DPangW_e'].Integral()
print "mu_a",  h['DPangW_mu'].Integral()
print "pi_a",  h['DPangW_pi'].Integral()
print "ka_a",  h['DPangW_ka'].Integral()
print "mix_a", h['DPangW_mix'].Integral()
print "oth_a", h['DPangW_oth'].Integral()
print "Reco Eff", h['DPangW'].Integral()
print "Vessel ", h['DPves'].Integral() 
print "e ", h['DP_e'].Integral()
print "mu ", h['DP_mu'].Integral()
print "tau ", h['DP_tau'].Integral()
print "pi ", h['DP_pi'].Integral()
print "ka ", h['DP_ka'].Integral()
print "oth ", h['DP_oth'].Integral()
print "mix ", h['DP_mix'].Integral()

o1 = tmp1+"_Ana_all.dat"
o2 = tmp1+"_Ana_e.dat"
o3 = tmp1+"_Ana_mu.dat"
o4 = tmp1+"_Ana_pi.dat"
o5 = tmp1+"_Ana_ka.dat"
o6 = tmp1+"_Ana_sum.dat"
o7 = tmp1+"_Ana_rate1.dat"
o8 = tmp1+"_Ana_rate2.dat"
o9 = tmp1+"_Ana_tau.dat" 
o0 = tmp1+"_Ana_mix.dat"
a=open(o1,'w+')
b=open(o2,'w+')
c=open(o3,'w+')
d=open(o4,'w+')
e=open(o5,'w+')
f=open(o6,'w+')
g=open(o7,'w+')
H=open(o8,'w+') 
k=open(o9,'w+')
l=open(o0,'w+')

DP_instance=darkphoton.DarkPhoton(float(mass_mc),float(eps))
other=0.
ves_other=0.
ang_other=0.
Nomv,Denv = 0., 0.
if float(h['DPvesW_e'].Integral())>0.:
    other+=float(h['DP_e'].Integral())
    ves_other+=float(h['DPves_e'].Integral())
    ang_other+=float(h['DPangW_e'].Integral())
    b.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, float(h['DP_e'].Integral())/float(h['DP'].Integral()), float(h['DPves_e'].Integral())/float(h['DP_e'].Integral()), float(h['DPangW_e'].Integral())/float(h['DPves_e'].Integral())))
    b.write('\n')
if float(h['DPvesW_mu'].Integral())!=0:
    other+=float(h['DP_mu'].Integral())
    ves_other+=float(h['DPves_mu'].Integral())
    ang_other+=float(h['DPangW_mu'].Integral())  
    c.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, float(h['DP_mu'].Integral())/float(h['DP'].Integral()), float(h['DPves_mu'].Integral())/float(h['DP_mu'].Integral()), float(h['DPangW_mu'].Integral())/float(h['DPves_mu'].Integral())))
    c.write('\n')
if float(h['DPvesW_tau'].Integral())!=0:
    other+=float(h['DP_tau'].Integral())
    ves_other+=float(h['DPves_tau'].Integral())
    ang_other+=float(h['DPangW_tau'].Integral())  
    k.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, float(h['DP_tau'].Integral())/float(h['DP'].Integral()),  float(h['DPves_tau'].Integral())/float(h['DP_tau'].Integral()), float(h['DPangW_tau'].Integral())/float(h['DPves_tau'].Integral())))
    k.write('\n')
if float(h['DPves_pi'].Integral())!=0:
    other+=float(h['DP_pi'].Integral())
    ves_other+=float(h['DPves_pi'].Integral())
    ang_other+=float(h['DPangW_pi'].Integral()) 
    Nomv+=float(h['DPves_pi'].Integral())
    Denv+=float(h['DP_pi'].Integral()) 
    d.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, float(h['DP_pi'].Integral())/float(h['DP'].Integral()), float(h['DPves_pi'].Integral())/float(h['DP_pi'].Integral()), float(h['DPangW_pi'].Integral())/float(h['DPves_pi'].Integral())))
    d.write('\n')
if float(h['DPvesW_ka'].Integral())!=0:
    other+=float(h['DP_ka'].Integral())
    ves_other+=float(h['DPves_ka'].Integral())
    ang_other+=float(h['DPangW_ka'].Integral()) 
    Nomv+=float(h['DPves_ka'].Integral())
    Denv+=float(h['DP_ka'].Integral())
    e.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, float(h['DP_ka'].Integral())/float(h['DP'].Integral()), float(h['DPves_ka'].Integral())/float(h['DP_ka'].Integral()), float(h['DPangW_ka'].Integral())/float(h['DPves_ka'].Integral())))
    e.write('\n')
if float(h['DPves_mix'].Integral())!=0: 
    other+=float(h['DP_mix'].Integral())
    ves_other+=float(h['DPves_mix'].Integral())
    ang_other+=float(h['DPangW_mix'].Integral()) 
    Nomv+=float(h['DPves_mix'].Integral())
    Denv+=float(h['DP_mix'].Integral())
    l.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, float(h['DP_mix'].Integral())/float(h['DP'].Integral()), float(h['DPves_mix'].Integral())/float(h['DP_mix'].Integral()), float(h['DPangW_mix'].Integral()/float(h['DPves_mix'].Integral()))))
    l.write('\n')
if float(ves_other)!=0.: 
    f.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, float(other)/float(h['DP'].Integral()), float(ves_other)/float(other), float(ang_other)/float(ves_other)))
    f.write('\n')
if float(h['DP'].Integral())!=0:
    a.write('%.4g %s %.8g %.8g %.8g %.8g %.8g' %(mass_mc, eps, nEvents, float(h['DP'].Integral()), float(other), float(h['DPves'].Integral()), float(h['DPangW'].Integral())))
    a.write('\n')

Nom   = 0.
Den   = 0.
NomL  = 0.
DenL  = 0.
 
if float(h['DP_e'].Integral())!=0:
    Nom+=float(h['DPang_e'].Integral())
    Den+=float(h['DP_e'].Integral())
    BR2=DP_instance.findBranchingRatio('A -> e- e+')
    NomL+=float(h['DPang_e'].Integral())
    DenL+=float(h['DP_e'].Integral())*BR2
if float(h['DP_mu'].Integral())!=0:
    Nom+=float(h['DPang_mu'].Integral())
    Den+=float(h['DP_mu'].Integral())
    BR1=DP_instance.findBranchingRatio('A -> mu- mu+')
    NomL+=float(h['DPang_mu'].Integral())
    DenL+=float(h['DP_mu'].Integral())*BR1
if float(h['DP_tau'].Integral())!=0:
    Nom+=float(h['DPang_tau'].Integral())
    Den+=float(h['DP_tau'].Integral())
    BR3=DP_instance.findBranchingRatio('A -> tau- tau+')
    NomL+=float(h['DPang_tau'].Integral())
    DenL+=float(h['DP_tau'].Integral())*BR3
if float(h['DP_pi'].Integral())!=0:
    Nom+=float(h['DPang_pi'].Integral())
    Den+=float(h['DP_pi'].Integral())
if float(h['DP_ka'].Integral())!=0:
    Nom+=float(h['DPang_ka'].Integral())
    Den+=float(h['DP_ka'].Integral())
if float(h['DP_mix'].Integral())!=0:
    Nom+=float(h['DPang_mix'].Integral())
    Den+=float(h['DP_mix'].Integral())
if float(h['DP'].Integral())!=0:
    RecLW=NomL/DenL*2.0e+20
    RecW=Nom/h['DP'].Integral()*2.0e+20
    g.write('%.4g %s %.8g' %(mass_mc, eps, RecW)) 
    g.write('\n')       
    H.write('%.4g %s %.8g' %(mass_mc, eps, RecLW)) 
    H.write('\n')

a.close()
b.close()
c.close()
d.close()
e.close()
f.close()
g.close()
H.close()
k.close()
l.close()
hfile =tmp1+"_ana.root" 
r.gROOT.cd()#gr.
ut.writeHists(h,hfile)
print hfile
print h['DP'].Integral() 
print RecW
print RecLW
