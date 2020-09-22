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
dpMom = '' 

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:p:m:e:A:g:f:", ["date=","production=","mass=","epsilon=","motherID=","geoFile=","final_dest="])
except getopt.GetoptError:
    print 'no file'
    sys.exit()
for o,a in opts:
    if o in ('-d',): date = a
    if o in ('-p',): pro = a
    if o in ('-m',): mass_mc = a
    if o in ('-e',): eps = a
    if o in ('-A',): dpMom = a
    if o in ('-g', '--geoFile',): geoFile = a
    if o in ('-f',): dest = a

if dpMom!='': tmp1 = "/eos/experiment/ship/data/DarkPhoton/PBC-June-3/"+date+"/reco/"+pro+"_"+dpMom+"_mass"+mass_mc+"_eps"+eps
if dpMom=='': tmp1 = "/eos/experiment/ship/data/DarkPhoton/PBC-June-3/"+date+"/reco/"+pro+"_mass"+mass_mc+"_eps"+eps
if pro=="pbrem1": tmp1 = "/eos/experiment/ship/data/DarkPhoton/PBC-June-3/"+date+"/reco/pbrem_mass"+mass_mc+"_eps"+eps
inputFile = tmp1+"_rec.root"
print inputFile
mass_mc=float(mass_mc)
eps=float(eps)
 

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

PDG = r.TDatabasePDG.Instance()

def dist2InnerWall(X,Y,Z):
    dist = 0
    node = sGeo.FindNode(X,Y,Z)
    #print node.GetName()
    if ShipGeo.tankDesign < 5:
        if not 'cave' in node.GetName(): return dist  # TP 
    else:
        if not 'decayVol' in node.GetName(): return dist
        #if not 'DecayVolume' in node.GetName(): return dist
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
    #print minDistance
    return minDistance

def isInFiducial(X,Y,Z):
    if Z > ShipGeo.TrackStation1.z : return False
    if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
    if dist2InnerWall(X,Y,Z)<5*u.cm: return False
    return True 

def findmum():#this function finds the mother of DP with weight,xs,momentum etc. USED for finding DP event
    for dp_ind,dp_tr in enumerate(sTree.MCTrack):
        if dp_tr.GetPdgCode()==9900015 or dp_tr.GetPdgCode()==4900023:
            mum_id=dp_tr.GetMotherId()
            #print mum_id
            dp_id=dp_ind
            #print dp_id
            if pro=='qcd' and dp_id==0: continue
            #print mum_id 
            mum_pdg=sTree.MCTrack[mum_id].GetPdgCode()
            #print mum_pdg
            if pro=='meson':
                xsw = dputil.getDPprodRate(mass_mc,eps,'meson',mum_pdg)
                if 'eta1' in dpMom and xsw!=0:
                    xsw1=xsw[1]
                    xsw=xsw[0]
            else: xsw = dputil.getDPprodRate(mass_mc,eps,pro,0) 
            #print "bu da farkli", xsw
            wg = sTree.MCTrack[dp_id].GetWeight()
            #print wg
            #print dp_id 
            dp_mom=r.TVector3(sTree.MCTrack[dp_id].GetPx(),sTree.MCTrack[dp_id].GetPy(),sTree.MCTrack[dp_id].GetPz())
            dp_mag=sTree.MCTrack[dp_id].GetP()
            break
        else:
            if 'eta1' in dpMom: xsw,xsw1,wg,dp_id,dp_mom,dp_mag,mum_id=0,0,0,0,0,0,0
            if not 'eta1' in dpMom: xsw,wg,dp_id,dp_mom,dp_mag,mum_id=0,0,0,0,0,0
    if 'eta1' in dpMom: return xsw,xsw1,wg,dp_id,dp_mom,dp_mag,mum_id
    if not 'eta1' in dpMom: return xsw,wg,dp_id,dp_mom,dp_mag,mum_id

def totCharge(pdg):
    CC=0
    PRT=PDG.GetParticle(pdg)
    if PRT.Charge()>0.:
        CC=+1
        return CC
    elif PRT.Charge()<0.:
        CC=-1
        return CC
    else: return CC

def findLepton(pdg):
    if abs(pdg)==11 or abs(pdg)==13 or abs(pdg)==15: return True
    else: return False

def checkLepMode(sTree, dp_id):
    PID=[]
    for mc,tr in enumerate(sTree.MCTrack):
        pid = tr.GetPdgCode()
        mom = tr.GetMotherId() 
        if mc>1:
            mom_pid=sTree.MCTrack[mom].GetPdgCode()
            if abs(pid)>9 and abs(pid)!=21 and findLepton(pid):
                if mom==dp_id:#leptons in all process and/or chargrons in meson process
                    PID.append(mc)
                elif(abs(mom_pid)<9 or abs(mom_pid)==21):#chargrons in qcd and pbrem
                    PID.append(mc)
                #else:
                    #if tr.GetProcID()==0. and findStable(pid): PID.append(mc)
            if tr.GetProcID()!=0.: return PID
    return PID

tmp1=tmp1.replace(date,dest)
tmp1=tmp1.replace("reco","ana")

if pro=='pbrem1':
    tmp1=tmp1.replace("pbrem","pbrem1")

Rfile=r.TFile(tmp1+"_miniShield.root",'recreate')
miniShield=r.TTree("miniShield","results of the mumu channel for mini-shield study")

Mu=r.std.vector(int)()
#Mc=r.std.vector(int)()
Ev=r.std.vector(int)()
Ves=r.std.vector(int)()
Wg=r.std.vector(float)()

miniShield.Branch('Mu',Mu)
#miniShield.Branch('Mc',Mc)
miniShield.Branch('Ves',Ves)
miniShield.Branch('Ev',Ev)
miniShield.Branch('Wg',Wg)

Vtx_x=r.std.vector(float)()
Vtx_y=r.std.vector(float)()
Vtx_z=r.std.vector(float)()
P_x=r.std.vector(float)()
P_y=r.std.vector(float)()
P_z=r.std.vector(float)()

Vtx_ves_x=r.std.vector(float)()
Vtx_ves_y=r.std.vector(float)()
Vtx_ves_z=r.std.vector(float)()
P_ves_x=r.std.vector(float)()
P_ves_y=r.std.vector(float)()
P_ves_z=r.std.vector(float)()

Vtx_W_x =r.std.vector(float)()
Vtx_W_y =r.std.vector(float)()
Vtx_W_z =r.std.vector(float)()
P_W_x   =r.std.vector(float)()
P_W_y   =r.std.vector(float)()
P_W_z   =r.std.vector(float)()
 
Vtx_W_ves_x =r.std.vector(float)()
Vtx_W_ves_y =r.std.vector(float)()
Vtx_W_ves_z =r.std.vector(float)()
P_W_ves_x   =r.std.vector(float)()
P_W_ves_y   =r.std.vector(float)()
P_W_ves_z   =r.std.vector(float)()

miniShield.Branch('Vtx_x',Vtx_x)
miniShield.Branch('Vtx_y',Vtx_y)
miniShield.Branch('Vtx_z',Vtx_z)
miniShield.Branch('P_x',P_x)
miniShield.Branch('P_y',P_y)
miniShield.Branch('P_z',P_z)

miniShield.Branch('Vtx_ves_x',Vtx_ves_x)
miniShield.Branch('Vtx_ves_y',Vtx_ves_y)
miniShield.Branch('Vtx_ves_z',Vtx_ves_z)
miniShield.Branch('P_ves_x',P_ves_x)
miniShield.Branch('P_ves_y',P_ves_y)
miniShield.Branch('P_ves_z',P_ves_z)

miniShield.Branch('Vtx_W_x',Vtx_W_x)
miniShield.Branch('Vtx_W_y',Vtx_W_y)
miniShield.Branch('Vtx_W_z',Vtx_W_z)
miniShield.Branch('P_W_x',P_W_x)
miniShield.Branch('P_W_y',P_W_y)
miniShield.Branch('P_W_z',P_W_z)

miniShield.Branch('Vtx_W_ves_x',Vtx_W_ves_x)
miniShield.Branch('Vtx_W_ves_y',Vtx_W_ves_y)
miniShield.Branch('Vtx_W_ves_z',Vtx_W_ves_z)
miniShield.Branch('P_W_ves_x',P_W_ves_x)
miniShield.Branch('P_W_ves_y',P_W_ves_y)
miniShield.Branch('P_W_ves_z',P_W_ves_z)

h={}

ut.bookHist(h,'Vx_W','',100,-10000.,10000.)
ut.bookHist(h,'Vy_W','',100,-10000.,10000.)
ut.bookHist(h,'Vz_W','',100,-9800.,4200.)
ut.bookHist(h,'Vx_Wv','',100,-10000.,10000.)
ut.bookHist(h,'Vy_Wv','',100,-10000.,10000.)
ut.bookHist(h,'Vz_Wv','',100,-9800.,4200.)
ut.bookHist(h,'Vx','',100,-10000.,10000.)
ut.bookHist(h,'Vy','',100,-10000.,10000.)
ut.bookHist(h,'Vz','',100,-9800.,4200.)
ut.bookHist(h,'Vx_v','',100,-10000.,10000.)
ut.bookHist(h,'Vy_v','',100,-10000.,10000.)
ut.bookHist(h,'Vz_v','',100,-9800.,4200.)

ut.bookHist(h,'Px_W','',100,0.,400.)
ut.bookHist(h,'Py_W','',100,0.,400.)
ut.bookHist(h,'Pz_W','',100,0.,400.)
ut.bookHist(h,'Px_Wv','',100,0.,400.)
ut.bookHist(h,'Py_Wv','',100,0.,400.)
ut.bookHist(h,'Pz_Wv','',100,0.,400.)
ut.bookHist(h,'Px','',100,0.,400.)
ut.bookHist(h,'Py','',100,0.,400.)
ut.bookHist(h,'Pz','',100,0.,400.)
ut.bookHist(h,'Px_v','',100,0.,400.)
ut.bookHist(h,'Py_v','',100,0.,400.)
ut.bookHist(h,'Pz_v','',100,0.,400.)

def myEventLoop(n):# Analysis is starting here
    Mu.clear()
    #Mc.clear()
    Ev.clear()
    Wg.clear()
    Ves.clear()

    Vtx_x.clear() 
    Vtx_y.clear()
    Vtx_z.clear()
    P_x.clear()
    P_y.clear()
    P_z.clear()
    Vtx_ves_x.clear() 
    Vtx_ves_y.clear()
    Vtx_ves_z.clear()
    P_ves_x.clear()
    P_ves_y.clear()
    P_ves_z.clear()
    Vtx_W_x.clear()
    Vtx_W_y.clear()
    Vtx_W_z.clear()
    P_W_x.clear()
    P_W_y.clear()
    P_W_z.clear()
    Vtx_W_ves_x.clear()
    Vtx_W_ves_y.clear()
    Vtx_W_ves_z.clear()
    P_W_ves_x.clear()
    P_W_ves_y.clear()
    P_W_ves_z.clear()

    rc=sTree.GetEntry(n)
    fm=findmum()
    if 'eta1' in dpMom:
        xsw=fm[0]
        xsw1=fm[1]
        wg=fm[2]
        dp_id=fm[3]
        dp_M=fm[4]
        dp_Mag=fm[5]
        mum=fm[6]
    if not 'eta1' in dpMom:
        xsw=fm[0]
        wg=fm[1]
        dp_id=fm[2]
        dp_M=fm[3]
        dp_Mag=fm[4]
        mum=fm[5]
    mu =0
    vessel = 0
    CM = 0
    vtxX, vtxY, vtxZ = [], [], []
    pX, pY, pZ = [], [], []
    if xsw==0 and wg==0 and dp_id==0: 
        #Dump(sTree.MCTrack)
        return 0
    dau=checkLepMode(sTree,dp_id) 
    for xxx in dau:
        if abs(sTree.MCTrack[xxx].GetPdgCode())==13:
            mu+=1
            CM+=totCharge(sTree.MCTrack[xxx].GetPdgCode())
            vtxX.append(sTree.MCTrack[xxx].GetStartX())
            vtxY.append(sTree.MCTrack[xxx].GetStartY())
            vtxZ.append(sTree.MCTrack[xxx].GetStartZ())
            pX.append(sTree.MCTrack[xxx].GetPx())
            pY.append(sTree.MCTrack[xxx].GetPy())
            pZ.append(sTree.MCTrack[xxx].GetPz())
            if isInFiducial(sTree.MCTrack[xxx].GetStartX(),sTree.MCTrack[xxx].GetStartY(),sTree.MCTrack[xxx].GetStartZ()):  vessel+=1

    if CM==0:
        Ev.push_back(n)
        Wg.push_back(wg)
        Mu.push_back(mu)
        #Mc.push_back(CM)
        for m in range(mu):
            Ves.push_back(vessel)
            Vtx_x.push_back(vtxX[m])
            Vtx_y.push_back(vtxY[m])
            Vtx_z.push_back(vtxZ[m])
            P_x.push_back(pX[m])
            P_y.push_back(pY[m])
            P_z.push_back(pZ[m])
            Vtx_W_x.push_back(vtxX[m]*wg)
            Vtx_W_y.push_back(vtxY[m]*wg)
            Vtx_W_z.push_back(vtxZ[m]*wg)
            P_W_x.push_back(pX[m]*wg)
            P_W_y.push_back(pY[m]*wg)
            P_W_z.push_back(pZ[m]*wg)
            if vessel:
                Vtx_ves_x.push_back(vtxX[m])
                Vtx_ves_y.push_back(vtxY[m])
                Vtx_ves_z.push_back(vtxZ[m])
                P_ves_x.push_back(pX[m])
                P_ves_y.push_back(pY[m])
                P_ves_z.push_back(pZ[m])
                Vtx_W_ves_x.push_back(vtxX[m]*wg)
                Vtx_W_ves_y.push_back(vtxY[m]*wg)
                Vtx_W_ves_z.push_back(vtxZ[m]*wg)
                P_W_ves_x.push_back(pX[m]*wg)
                P_W_ves_y.push_back(pY[m]*wg)
                P_W_ves_z.push_back(pZ[m]*wg)
        if mu==2:
            h['Vx'].Fill(vtxX[0])
            h['Vy'].Fill(vtxY[0])
            h['Vz'].Fill(vtxZ[0])
            h['Px'].Fill(pX[0])
            h['Py'].Fill(pY[0])
            h['Pz'].Fill(pZ[0])
            h['Vx'].Fill(vtxX[1])
            h['Vy'].Fill(vtxY[1])
            h['Vz'].Fill(vtxZ[1])
            h['Px'].Fill(pX[1])
            h['Py'].Fill(pY[1])
            h['Pz'].Fill(pZ[1])
            h['Vx_W'].Fill(vtxX[0],wg) 
            h['Vy_W'].Fill(vtxY[0],wg)
            h['Vz_W'].Fill(vtxZ[0],wg)
            h['Px_W'].Fill(pX[0],wg)
            h['Py_W'].Fill(pY[0],wg)
            h['Pz_W'].Fill(pZ[0],wg)
            h['Vx_W'].Fill(vtxX[1],wg)
            h['Vy_W'].Fill(vtxY[1],wg)
            h['Vz_W'].Fill(vtxZ[1],wg)
            h['Px_W'].Fill(pX[1],wg)
            h['Py_W'].Fill(pY[1],wg)
            h['Pz_W'].Fill(pZ[1],wg)
            if vessel==2:
                h['Vx_v'].Fill(vtxX[0])
                h['Vy_v'].Fill(vtxY[0])
                h['Vz_v'].Fill(vtxZ[0])
                h['Px_v'].Fill(pX[0])
                h['Py_v'].Fill(pY[0])
                h['Pz_v'].Fill(pZ[0])
                h['Vx_v'].Fill(vtxX[1])
                h['Vy_v'].Fill(vtxY[1])
                h['Vz_v'].Fill(vtxZ[1])
                h['Px_v'].Fill(pX[1])
                h['Py_v'].Fill(pY[1])
                h['Pz_v'].Fill(pZ[1])
                h['Vx_Wv'].Fill(vtxX[0],wg)
                h['Vy_Wv'].Fill(vtxY[0],wg)
                h['Vz_Wv'].Fill(vtxZ[0],wg)
                h['Px_Wv'].Fill(pX[0],wg)
                h['Py_Wv'].Fill(pY[0],wg)
                h['Pz_Wv'].Fill(pZ[0],wg)
                h['Vx_Wv'].Fill(vtxX[1],wg) 
                h['Vy_Wv'].Fill(vtxY[1],wg)
                h['Vz_Wv'].Fill(vtxZ[1],wg)
                h['Px_Wv'].Fill(pX[1],wg)
                h['Py_Wv'].Fill(pY[1],wg)
                h['Pz_Wv'].Fill(pZ[1],wg)
        miniShield.Fill()
        if vessel==1 and mu==2: print vessel, mu, n, CM
    #if mu>2: Dump(sTree.MCTrack)

nEvents =sTree.GetEntries()
for n in range(nEvents):
    myEventLoop(n)

Rfile.Write()
Rfile.Close()
r.gROOT.cd()
ut.writeHists(h,tmp1+"_results.root")
