import ROOT
import os,sys,getopt
from decorators import *
import rootUtils as ut
import shipunit as u
from rootpyPickler import Unpickler
import shipRoot_conf
#import alpProductionRates as alputil
import math as m
import numpy as np
from array import array
import subprocess
shipRoot_conf.configure()
sim=0

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:m:G:g:f:s:", ["date=","mass=","coupling=","geoFile=","final_dest=","sim="])

except getopt.GetoptError:
    print 'no file'
    sys.exit()

for o,a in opts:
    if o in ('-d',): date = a
    if o in ('-m',): mass_mc = a
    if o in ('-G',): gax = a
    if o in ('-g', '--geoFile',): geoFile = a
    if o in ('-f',): dest = a
    if o in ('-s',): sim = int(a)

tmp1 = "/eos/experiment/ship/data/ALP/"+date+"/runs/ALPACA_mass"+mass_mc+"_g"+gax
gax=gax.replace(".000000","")
pipe = subprocess.Popen("ls /eos/experiment/ship/data/ALP/"+date+"/ntuples/alp_m"+mass_mc+"_g"+gax+"*.root", shell=True, stdout=subprocess.PIPE).stdout
xs=pipe.read()
xs=xs.replace("/eos/experiment/ship/data/ALP/"+date+"/ntuples/alp_m"+mass_mc+"_g"+gax+"_xs","")
xs=xs.replace(".root","")
print xs
#inputFile = tmp1+".root"
mass_mc=float(mass_mc)
gax=float(gax)

eosship =  ROOT.gSystem.Getenv("EOSSHIP")
eospath = eosship+geoFile
#fgeo = ROOT.TFile.Open(eospath)
fgeo = ROOT.TFile.Open(geoFile)
#sGeo = ROOT.gGeoManager

upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')
dy = ShipGeo.Yheight/u.m
MeasCut=25

# -----Create geometry----------------------------------------------
import shipDet_conf
run = ROOT.FairRunSim()
run.SetName("TGeant4")  # Transport engine
run.SetOutputFile(ROOT.TMemFile('output', 'recreate'))  # Output file
run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for the mag field
# -----Create geometry----------------------------------------------
modules = shipDet_conf.configure(run,ShipGeo)

import geomGeant4
if hasattr(ShipGeo.Bfield,"fieldMap"):
    fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True, withVirtualMC = False)
else:
    print("no fieldmap given, geofile too old, not anymore support")
    exit(-1)

sGeo   = fgeo.FAIRGeom
bfield = ROOT.genfit.FairShipFields()
bfield.setField(fieldMaker.getGlobalField())
fM = ROOT.genfit.FieldManager.getInstance()
fM.init(bfield)
PDG = ROOT.TDatabasePDG.Instance()

import TrackExtrapolateTool

def dist2InnerWall(X,Y,Z):
    dist = 0
    # return distance to inner wall perpendicular to z-axis, if outside decayVolume return 0.
    node = sGeo.FindNode(X,Y,Z)
    if ShipGeo.tankDesign < 5:
        if not 'cave' in node.GetName(): return dist  # TP 
    else:
        if not 'DecayVacuum' in node.GetName(): return dist
    start = array('d',[X,Y,Z])
    nsteps = 8
    dalpha = 2*ROOT.TMath.Pi()/nsteps
    rsq = X**2+Y**2
    minDistance = 100 *u.m
    for n in range(nsteps):
        alpha = n * dalpha
        sdir  = array('d',[ROOT.TMath.Sin(alpha),ROOT.TMath.Cos(alpha),0.])
        node = sGeo.InitTrack(start, sdir)
        nxt = sGeo.FindNextBoundary()
        if ShipGeo.tankDesign < 5 and nxt.GetName().find('I')<0: return 0    
        distance = sGeo.GetStep()
        if distance < minDistance  : minDistance = distance
    return minDistance

def checkFiducialVolume(sTree,xxx,dy):
# extrapolate track to middle of magnet and check if in decay volume
    inside = True
    rc,pos,mom = False,None,None
    parallelToZ = ROOT.TVector3(0., 0., 1.)
    NewPosition = ROOT.TVector3(0., 0., ShipGeo.SplitCal.ZStart) 
    rep    = ROOT.genfit.RKTrackRep(13*cmp(sTree.MCTrack[xxx].GetPdgCode(),0) )
    state  = ROOT.genfit.StateOnPlane(rep)
    print state
    vtx0=ROOT.TVector3()
    mom0=ROOT.TVector3()
    sTree.MCTrack[xxx].GetStartVertex(vtx0)
    sTree.MCTrack[xxx].GetMomentum(mom0)
    pos,mom = vtx0, mom0
    rep.setPosMom(state,pos,mom)
    try:    
        rep.extrapolateToPlane(state, NewPosition, parallelToZ )
        pos,mom = state.getPos(),state.getMom()
        print pos.X(), pos.Y(), pos.Z()
        rc = True 
    except: 
        # print 'error with extrapolation: z=',z/u.m,'m',pos.X(),pos.Y(),pos.Z(),mom.X(),mom.Y(),mom.Z()
        pass 
    if not rc:
         # use linear extrapolation
        px,py,pz  = mom.X(),mom.Y(),mom.Z()
        lam = (ShipGeo.SplitCal.ZStart-pos.Z())/pz
        pos = ROOT.TVector3( pos.X()+lam*px, pos.Y()+lam*py, ShipGeo.SplitCal.ZStart )
    if not rc: return False
    if not dist2InnerWall(pos.X(),pos.Y(),pos.Z())>0: return False
    return inside

def isInFiducial(X,Y,Z):
    if Z > ShipGeo.TrackStation1.z : return False
    if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
    if dist2InnerWall(X,Y,Z)<5*u.cm: return False
    return True 

def findMom():#this function finds the mother of ALP with weight,xs,momentum etc. USED for finding ALP event
    for ind, tr in enumerate(sTree.MCTrack):
        if tr.GetPdgCode()==9900015:
            alp_ind = ind
            #xsw = dputil.getALPprodRate(mass_mc,eps)
            xsw=float(xs)*10**-9
            wg = sTree.MCTrack[alp_ind].GetWeight()
            alp_P=ROOT.TVector3(sTree.MCTrack[alp_ind].GetPx(),sTree.MCTrack[alp_ind].GetPy(),sTree.MCTrack[alp_ind].GetPz())
            alp_Mag=sTree.MCTrack[alp_ind].GetP()
            #print sTree.MCTrack[alp_ind].GetStartZ()
            break
        else: xsw, wg, alp_ind, alp_P, alp_Mag = 0, 0, 0, 0, 0
    return  xsw, wg, alp_ind, alp_P, alp_Mag

def find_photon(pdg):
    if pdg==22: return True
    else: return False

def findDau(sTree, alp_ind):
    PID=[]
    for mc,tr in enumerate(sTree.MCTrack):
        pid = tr.GetPdgCode()
        mom = tr.GetMotherId()
        if mom==alp_ind and find_photon(pid):
            PID.append(mc)
        if len(PID)==2:return PID

def firstProd(sTree,dauPID):
    for mc,tr in enumerate(sTree.MCTrack):
        #if tr.GetMotherId()==dauPID and abs(tr.GetPdgCode())==11:
        if tr.GetMotherId()==dauPID:
            vtx=ROOT.TVector3()
            tr.GetStartVertex(vtx)
            #if abs(tr.GetPdgCode())==11 and vtx.Z()<=ShipGeo.SplitCal.SplitCalThickness+ShipGeo.SplitCal.ZStart and  vtx.Z()>=ShipGeo.SplitCal.ZStart and abs(vtx.X())<=ShipGeo.SplitCal.XMax and abs(vtx.Y())<=ShipGeo.SplitCal.YMax: return True
            if abs(tr.GetPdgCode())==11 and vtx.Z()<=ShipGeo.SplitCal.SplitCalThickness+ShipGeo.SplitCal.ZStart and  vtx.Z()>=ShipGeo.TrackStation1.z and abs(vtx.X())<=ShipGeo.SplitCal.XMax and abs(vtx.Y())<=ShipGeo.SplitCal.YMax: return True
    return False

def MCextrapolate(vtx,mom):
    lam=(ShipGeo.SplitCal.ZStart-vtx.Z())/mom.Z()
    if abs(vtx.X()+lam*mom.X())<=ShipGeo.SplitCal.XMax and abs(vtx.Y()+lam*mom.Y())<=ShipGeo.SplitCal.YMax: return True

def findSignal(sTree,xxx):
    PID=[]
    for mc,tr in enumerate(sTree.MCTrack):
        if tr.GetMotherId()==xxx and abs(tr.GetPdgCode())==11: PID.append(mc)
    return PID
def inDet(vtx):
    if vtx.Z()<=ShipGeo.SplitCal.SplitCalThickness+ShipGeo.SplitCal.ZStart and  vtx.Z()>=ShipGeo.SplitCal.ZStart and abs(vtx.X())<=ShipGeo.SplitCal.XMax and abs(vtx.Y())<=ShipGeo.SplitCal.YMax: return True
    else: return False

def IPtoTarget(vtx, mom, Avtx):
    p = mom.P()
    delta = 0.
    for i in range(3):
        delta += (Avtx(i) - vtx(i)) * mom(i)/p
    ip = 0.
    for i in range(3):
        ip += (Avtx(i) - vtx(i) - delta*mom(i)/p)**2.
    return ROOT.TMath.Sqrt(ip)
h={}

ut.bookHist(h,'ALPangWxs','invariant Mass with Weights (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALP','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPW','invariant Mass with Weights (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPpur','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPpurW','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPves','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPvesW','invariant Mass with Weights (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPang','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPangW','invariant Mass with Weights (GeV)',100,0.,mass_mc+5.)

ut.bookHist(h,'ALP_oth','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPW_oth','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPpur_oth','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPpurW_oth','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPves_oth','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPvesW_oth','invariant Mass with Weights (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPang_oth','invariant Mass (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'ALPangW_oth','invariant Mass with Weights (GeV)',100,0.,mass_mc+5.)
ut.bookHist(h,'IP','',100,0.,1000.)

def myEventLoop(n):# Analysis is starting here
    DAU=0
    VES=0
    REC=0
    rc=sTree.GetEntry(n) 
    fm=findMom()
    xsw=fm[0]
    wg=fm[1]
    alp_ind=fm[2]
    alp_P=fm[3]
    alp_Mag=fm[4]
    if xsw==0 and wg==0 and alp_ind==0:
        h['ALP_oth'].Fill(mass_mc) 
        h['ALPW_oth'].Fill(mass_mc,wg)
        return 0
    #Dump(sTree.MCTrack)
    h['ALP'].Fill(mass_mc) 
    h['ALPW'].Fill(mass_mc,wg)
    dau=findDau(sTree,alp_ind)
    for xxx in dau:
        if sTree.MCTrack[xxx].GetPdgCode()==22:
            DAU+=1
            vtx=ROOT.TVector3(sTree.MCTrack[xxx].GetStartX(), sTree.MCTrack[xxx].GetStartY(), sTree.MCTrack[xxx].GetStartZ())
            mom=ROOT.TVector3(sTree.MCTrack[xxx].GetPx(), sTree.MCTrack[xxx].GetPy(), sTree.MCTrack[xxx].GetPz())
            Mom=ROOT.TLorentzVector()
            sTree.MCTrack[xxx].Get4Momentum(Mom)
            Vtx=ROOT.TLorentzVector(sTree.MCTrack[xxx].GetStartX(), sTree.MCTrack[xxx].GetStartY(), sTree.MCTrack[xxx].GetStartZ(), sTree.MCTrack[xxx].GetStartT())
            Avtx=ROOT.TLorentzVector(sTree.MCTrack[0].GetStartX(), sTree.MCTrack[0].GetStartY(), sTree.MCTrack[0].GetStartZ(), sTree.MCTrack[0].GetStartT())
            if isInFiducial(vtx.X(),vtx.Y(),vtx.Z()):
                VES+=1
                IP=IPtoTarget(Vtx, Mom,Avtx)
                h['IP'].Fill(IP)
                if firstProd(sTree,xxx) and sTree.MCTrack[0].GetStartT()-sTree.MCTrack[xxx].GetStartT()<=1. and MCextrapolate(vtx,mom) and mom.Mag()>=1.:
                    REC+=1
    """for tr, rec in enumerate(sTree.Reco_SplitcalClusters):
        e=+rec.GetEnergy()"""

    if DAU>1:
        h['ALPpur'].Fill(mass_mc) 
        h['ALPpurW'].Fill(mass_mc,wg)
        if VES>1:
            h['ALPves'].Fill(mass_mc) 
            h['ALPvesW'].Fill(mass_mc,wg)
            if REC>1:
                h['ALPang'].Fill(mass_mc) 
                h['ALPangW'].Fill(mass_mc,wg)
                h['ALPangWxs'].Fill(mass_mc,wg*xsw)
            elif REC<2:
                h['ALPang_oth'].Fill(mass_mc) 
                h['ALPangW_oth'].Fill(mass_mc,wg)
            #else: print("WRONG REC")
        elif VES<2:
            h['ALPves_oth'].Fill(mass_mc) 
            h['ALPvesW_oth'].Fill(mass_mc,wg)
        else: print("WRONG VES")
    elif DAU<2:
        h['ALPpur_oth'].Fill(mass_mc) 
        h['ALPpurW_oth'].Fill(mass_mc,wg)
    else: print("WRONG")

if not sim:
    for i in range(1,6):
        inputFile = tmp1+"_"+str(i)+".root"
        print inputFile
        eospath = eosship+inputFile
        try:
            f = ROOT.TFile.Open(eospath)
            sTree=f.cbmsim
            nEvents =sTree.GetEntries()
            #nEvents=100
            for n in range(nEvents):
                myEventLoop(n)
        except: continue

if sim:
    inputFile = tmp1+".root"
    print inputFile
    eospath = eosship+inputFile
    f = ROOT.TFile.Open(eospath)
    sTree=f.cbmsim
    nEvents =sTree.GetEntries()
    #nEvents=100
    for n in range(nEvents):
        myEventLoop(n)


print(h['ALP'].Integral(), h['ALPW'].Integral(), h['ALPpur'].Integral(), h['ALPpurW'].Integral(), h['ALPves'].Integral(), h['ALPvesW'].Integral(), h['ALPang'].Integral(), h['ALPangW'].Integral())
print(h['ALP_oth'].Integral(), h['ALPW_oth'].Integral(), h['ALPpur_oth'].Integral(), h['ALPpurW_oth'].Integral(), h['ALPves_oth'].Integral(), h['ALPvesW_oth'].Integral(), h['ALPang_oth'].Integral(), h['ALPangW_oth'].Integral())

tmp2=tmp1.replace(date,dest)
#tmp2=tmp2.replace("hadd","ana")
#tmp1=tmp1.replace("hadd","ana/dat")
tmp2=tmp2.replace("runs","ana")
tmp1=tmp1.replace("runs","ana/dat")
tmp1=tmp1.replace(date,dest)

o1  = tmp1+"_ratio.dat"
o2  = tmp1+"_sum.dat"
o3  = tmp1+"_other.dat"
o4  = tmp1+"_rate.dat"

a=open(o1,'w+')
b=open(o2,'w+')
c=open(o3,'w+')
d=open(o4,'w+')

if float(h['ALP'].Integral())!=0.0:
    b.write('{} {} {} {} {} {} {}'.format(mass_mc, gax, nEvents, h['ALP'].Integral(), h['ALPpur'].Integral(), h['ALPves'].Integral(), h['ALPang'].Integral()))
    b.write('\n')
    if float(h['ALPves'].Integral())!=0.0:
        a.write('{} {} {} {} {} {}'.format(mass_mc, gax, h['ALP'].Integral()/nEvents, h['ALPpur'].Integral()/ h['ALP'].Integral(), h['ALPvesW'].Integral()/h['ALPpur'].Integral(), h['ALPangW'].Integral()/h['ALPvesW'].Integral()))
        a.write('\n')
    if float(h['ALPvesW'].Integral())==0.0:
        a.write('{} {} {} {} {} 0.'.format(mass_mc, gax, h['ALP'].Integral()/nEvents, h['ALPpur'].Integral()/ h['ALP'].Integral(), h['ALPvesW'].Integral()/h['ALPpur'].Integral()))
        a.write('\n')
    RecW=h['ALPangWxs'].Integral()/h['ALP'].Integral()*2.0e+20#weighted Selected/weighted Vessel
    d.write('{} {} {}'.format(mass_mc, gax, RecW))
    d.write('\n')

if float(h['ALP_oth'].Integral())!=0.0:
    if float(h['ALPves_oth'].Integral())!=0.0:
        c.write('{} {} {} {} {} {}'.format(mass_mc, gax, h['ALP_oth'].Integral()/nEvents, h['ALPpur_oth'].Integral()/ h['ALP_oth'].Integral(), h['ALPvesW_oth'].Integral()/h['ALPpur_oth'].Integral(), h['ALPangW_oth'].Integral()/h['ALPvesW_oth'].Integral()))
        c.write('\n')
    if float(h['ALPvesW_oth'].Integral())==0.0:
        c.write('{} {} {} {} {} 0.'.format(mass_mc, gax, h['ALP_oth'].Integral()/nEvents, h['ALPpur_oth'].Integral()/ h['ALP_oth'].Integral(), h['ALPvesW_oth'].Integral()/h['ALPpur_oth'].Integral()))
        c.write('\n')

print mass_mc, gax, RecW

a.close()
b.close()
c.close()
d.close()

tmp1=tmp1.replace("dat/","")
hfile =tmp2+"_ana.root" 
ROOT.gROOT.cd()
ut.writeHists(h,hfile)

