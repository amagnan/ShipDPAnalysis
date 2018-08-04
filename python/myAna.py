#What are left BR vs MAss, Eff vs Mass, output.dat for sensitivityPlot.C
import ROOT,os,sys,getopt
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

PDG = ROOT.TDatabasePDG.Instance()

try:
	opts, args = getopt.getopt(sys.argv[1:], "f:g:", ["nEvents=","geoFile="])

except getopt.GetoptError:
	print 'no file'
	sys.exit()

for o,a in opts:
	if o in ('-f',): inputFile = a
	if o in ('-g', '--geoFile',): geoFile = a

tmp=inputFile.replace('/eos/experiment/ship/data/DarkPhoton/PBC-June-3/rec/','')
tmp1=tmp.replace('_rec.root','')
tmp2=tmp1.replace('mass','')
tmp3=tmp2.replace('eps','')		
out=tmp3.split('_') 
pro=out[0]
mass_mc=float(out[1])
eps=float(out[2])


eosship =  ROOT.gSystem.Getenv("EOSSHIP")
eospath = eosship+inputFile
f = ROOT.TFile.Open(eospath)
sTree=f.cbmsim
eospath = eosship+geoFile
fgeo = ROOT.TFile.Open(eospath)

sGeo = ROOT.gGeoManager
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')

dy = ShipGeo.Yheight/u.m

import shipDet_conf

run = ROOT.FairRunSim()
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

geoMat =  ROOT.genfit.TGeoMaterialInterface()
ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
bfield = ROOT.genfit.FairShipFields()
fM = ROOT.genfit.FieldManager.getInstance()
fM.init(bfield)


h={}

ut.bookHist(h,'DP','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_tau','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DP_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)


ut.bookHist(h,'DPang','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_tau','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPang_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPangW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)


ut.bookHist(h,'DPves','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_e','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_e','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_mu','invariant Mass (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_mu','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_pi','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_pi','invariant Mass (GeV) with Weight',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_ka','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_ka','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_tau','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_tau','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPves_oth','invariant Mass',50,mass_mc-1.,mass_mc+1.)
ut.bookHist(h,'DPvesW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)

ut.bookHist(h,'eff_P0','',50,0.,250.)
ut.bookHist(h,'eff_P1','',50,0.,250.)
ut.bookHist(h,'eff_P2','',50,0.,250.)
ut.bookHist(h,'eff_Pt0','',50,0.,8.)
ut.bookHist(h,'eff_Pt1','',50,0.,8.)
ut.bookHist(h,'eff_Pt2','',50,0.,8.)
ut.bookHist(h,'eff_Eta0','',50,2.,8.)
ut.bookHist(h,'eff_Eta1','',50,2.,8.)
ut.bookHist(h,'eff_Eta2','',50,2.,8.)

def dist2InnerWall(X,Y,Z):
	dist = 0
	node = sGeo.FindNode(X,Y,Z)
	if ShipGeo.tankDesign < 5:
		if not 'cave' in node.GetName(): return dist  # TP 
	else:
		if not 'decayVol' in node.GetName(): return dist
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


import TrackExtrapolateTool

top = ROOT.gGeoManager.GetTopVolume()
ecal	= None
muon	= None
hcal	= None

if top.GetNode('Hcal_1'):
	hcal = top.GetNode('Hcal_1')
	z_hcal = hcal.GetMatrix().GetTranslation()[2]

if top.GetNode('MuonDetector_1'):
	muon = top.GetNode('MuonDetector_1')
	z_muon = muon.GetMatrix().GetTranslation()[2]

if top.GetNode('Ecal_1'):
	ecal = top.GetNode('Ecal_1')
	z_ecal = ecal.GetMatrix().GetTranslation()[2]


def isInFiducial(X,Y,Z):
	if Z > ShipGeo.TrackStation1.z : return False
	if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
	if dist2InnerWall(X,Y,Z)<5*u.cm: return False
	return True 


def checkFiducialVolume(fT,tkey,dy):
	inside = True
	fT = sTree.FitTracks[tkey]
	rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,ShipGeo.Bfield.z)#
	if not rc: return False
	return inside

def ang_muon(fT,dy):
	inside = True
	rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,z_muon)#
	if not rc: return False
	return inside


def ang_ecal(fT,dy):
	inside = True
	rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,z_ecal)#
	if not rc: return False
	return inside


def ang_hcal(fT,dy):
	inside = True
	rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,z_hcal)#
	if not rc: return False
	return inside


def checkDPorigin(sTree):
	flag= True
	dpkey=-1
	for n in range(sTree.MCTrack.GetEntries()):
		mo = sTree.MCTrack[n].GetMotherId()
		if mo<1: continue
		if abs(sTree.MCTrack[mo].GetPdgCode())==4900023 or abs(sTree.MCTrack[mo].GetPdgCode())==990015:
			dpkey=n
			break
		if sTree.MCTrack[mo].GetPdgCode()!= 4900023 and sTree.MCTrack[mo].GetPdgCode()!= 9900015 and abs(sTree.MCTrack[mo].GetPdgCode())!=11 and abs(sTree.MCTrack[mo].GetPdgCode())!=13 and (abs(mo.GetPdgCode())==4 or abs(sTree.MCTrack[mo].GetPdgCode())== 5 or abs(sTree.MCTrack[mo].GetPdgCode())==3 or abs(sTree.MCTrack[mo].GetPdgCode())== 2 or abs(sTree.MCTrack[mo].GetPdgCode())== 1) and (sTree.MCTrack[mo.GetMotherId()]==4900023 or sTree.MCTrack[mo.GetMotherId()]==9900015):
			dpkey=n
			break
	if dpkey<0:
		print 'no DP'
	else:
		dpVtx=sTree.MCTrack[dpkey]
		X,Y,Z = dpVtx.GetStartX(),dpVtx.GetStartY(),dpVtx.GetStartZ()
		if isInFiducial(X,Y,Z): flag = True
	return flag


def match2DP(p):
	matched=False
	dpKey=[]
	dpDau=[]
	artik=[]
	for t in [p.GetDaughter(0),p.GetDaughter(1)]:
		mcp=sTree.fitTrack2MC[t]
		#num=0
		while mcp> 0:
			mo=sTree.MCTrack[mcp]
			artik.append(abs(mo.GetPdgCode()))
			#num+=1
			if mo.GetPdgCode()== 4900023 or mo.GetPdgCode()==9900015:
				dpKey.append(mcp)
				dpDau.append(len(artik))
				break
			mcp=mo.GetMotherId()
	#if len(dpKey)==2 and artik[dpDau[0]-3]==artik[dpDau[1]-3]:
	if len(dpKey)==2 and len(artik)<7 and artik[dpDau[0]-3]==artik[dpDau[1]-3]:
		if dpKey[0]==dpKey[1]: matched=True
	if len(dpKey)==2 and artik[dpDau[0]-2]==artik[dpDau[1]-2] and artik[dpDau[0]-2]==11:
		if dpKey[0]==dpKey[1]: matched=True
	return matched

def myEventLoop(n,xsw):
	a_f, e_a, mu_a, pi_a, ka_a=0, 0, 0, 0, 0
	a_v, e_v, mu_v, pi_v, ka_v=0, 0, 0, 0, 0
	dp_pair,e_pair, mu_pair, pi_pair, ka_pair, neu_pair, oth_pair = 0, 0, 0, 0, 0, 0, 0
	rc=sTree.GetEntry(n)
	wg = sTree.MCTrack[1].GetWeight()
	#print n
	if not checkDPorigin(sTree): return#MCTRUE level; check vessel and mother
	#print n
	for f,fit in enumerate(sTree.FitTracks):
		#print f
		fitStatus = fit.getFitStatus()
		if fitStatus.isFitConverged():
			xx = fit.getFittedState()
			mc = sTree.MCTrack[sTree.fitTrack2MC[f]]
			vtx=ROOT.TVector3(mc.GetStartX(), mc.GetStartY(), mc.GetStartZ())
			mom=ROOT.TVector3(mc.GetPx(), mc.GetPy(), mc.GetPz())
			#print mc.GetPdgCode()
			if isInFiducial(vtx.X(),vtx.Y(),vtx.Z()):#IN checkDPorigin, it exist, but we  want to find right daughters
				a_v+=1
				if abs(mc.GetPdgCode())==11:
					e_v+=1
					if ang_ecal(fit,dy): e_a+=1
				if abs(mc.GetPdgCode())==15:
					tau_v+=1
					if ang_ecal(fit,dy): tau_a+=1
				if abs(mc.GetPdgCode())==13:
					mu_v+=1
					if ang_muon(fit,dy): mu_a+=1
				if abs(mc.GetPdgCode())==211:
					pi_v+=1
					if ang_ecal(fit,dy): pi_a+=1
				if abs(mc.GetPdgCode())==321:
					ka_v+=1
					if ang_ecal(fit,dy): ka_a+=1
				if checkFiducialVolume(sTree,f,dy): a_f+=1

	if a_f==2:

		dp_pair=1
		h['eff_P0'].Fill(sTree.MCTrack[1].GetP(),wg*xsw)			
		h['eff_Pt0'].Fill(sTree.MCTrack[1].GetPt(),wg*xsw)
		h['eff_Eta0'].Fill(sTree.MCTrack[1].GetRapidity(),wg*xsw)
		h['DPang'].Fill(mass_mc)
		h['DPangW'].Fill(mass_mc,wg*xsw)
		if e_a==2:
			e_pair=1
			h['DPang_e'].Fill(mass_mc)
			h['DPangW_e'].Fill(mass_mc,wg*xsw)
		if mu_a==2:
			mu_pair=1
			h['DPang_mu'].Fill(mass_mc)
			h['DPangW_mu'].Fill(mass_mc,wg*xsw)
		if pi_a==2:
			pi_pair=1
			h['DPang_pi'].Fill(mass_mc)
			h['DPangW_pi'].Fill(mass_mc,wg*xsw)
		if ka_a==2:
			ka_pair=1
			h['DPang_ka'].Fill(mass_mc)
			h['DPangW_ka'].Fill(mass_mc,wg*xsw)
	if (e_a==1 or mu_a==1) and (pi_a==1 or ka_a==1):
		oth_pair=1
		h['DPang_oth'].Fill(mass_mc)
		h['DPangW_oth'].Fill(mass_mc,wg*xsw)

	#DENOMINATOR OF Angular Acceptance
	if a_v==2:
		h['DPves'].Fill(mass_mc)
		h['DPvesW'].Fill(mass_mc,wg*xsw)
		if e_v==2:  
			h['DPves_e'].Fill(mass_mc)
			h['DPvesW_e'].Fill(mass_mc,wg*xsw)
		if mu_v==2:              
			h['DPves_mu'].Fill(mass_mc)
			h['DPvesW_mu'].Fill(mass_mc,wg*xsw)
		if pi_v==2:              
			h['DPves_pi'].Fill(mass_mc)
			h['DPvesW_pi'].Fill(mass_mc,wg*xsw)
		if ka_v==2:
			h['DPves_ka'].Fill(mass_mc)		
			h['DPvesW_ka'].Fill(mass_mc,wg*xsw)
	if (e_v==1 or mu_v==1) and (pi_v==1 or ka_v==1):
		h['DPves_oth'].Fill(mass_mc)
		h['DPvesW_oth'].Fill(mass_mc,wg*xsw)

	for r,DP in enumerate(sTree.Particles):
		t1,t2=DP.GetDaughter(0), DP.GetDaughter(1)#index of daughters
		if not checkFiducialVolume(sTree,t1,dy) or not checkFiducialVolume(sTree,t2,dy) : continue#ANGULAR ACCEPTANCE IN GENERAL?
		checkMeasurements = True
		for tr in [t1,t2]:
			fit=sTree.FitTracks[tr]
			fitStatus = fit.getFitStatus()
			if not fitStatus.isFitConverged(): continue
			xx = fit.getFittedState()
			nmeas = fitStatus.getNdf()
			if nmeas < 25: checkMeasurements = False
			if not checkMeasurements: continue
		DPpos=ROOT.TLorentzVector()
		DP.ProductionVertex(DPpos)
		DPmom=ROOT.TLorentzVector()
		DP.Momentum(DPmom)
		if not isInFiducial(DPpos.X(),DPpos.Y(),DPpos.Z()): continue#AFTER checkDPorigin IT MAY NOT BE NECCESARRY
		if not match2DP(DP): continue#FINAL RECO SELECTION
		mass = DPmom.M()
		#RATES AND GENERAL NOMINATOR
		#if r>0: print n,r
		h['eff_P1'].Fill(DP.P(),wg*xsw)
		h['eff_Pt1'].Fill(DP.Pt(),wg*xsw)
		h['eff_Eta1'].Fill(DP.Eta(),wg*xsw)
		if dp_pair==1:
			h['DP'].Fill(mass)
			h['DPW'].Fill(mass,wg*xsw)
		if e_pair==1:
			h['DP_e'].Fill(mass)
			h['DPW_e'].Fill(mass,wg*xsw)
		if mu_pair==1:
			h['DP_mu'].Fill(mass)
			h['DPW_mu'].Fill(mass,wg*xsw)
		if pi_pair==1:
			#if n>4000: Dump(sTree.MCTrack)
			h['DP_pi'].Fill(mass)
			h['DPW_pi'].Fill(mass,wg*xsw)
		if ka_pair==1:
			h['DP_ka'].Fill(mass)
			h['DPW_ka'].Fill(mass,wg*xsw)
		if oth_pair==1:
			h['DP_oth'].Fill(mass)
			h['DPW_oth'].Fill(mass,wg*xsw)


nEvents =sTree.GetEntries()
print nEvents
for n in range(nEvents):
	if pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,meson,9900015)#
	else: xsw = dputil.getDPprodRate(mass_mc,eps,pro,0)#i am taking exact mass?
	myEventLoop(n,xsw)
	sTree.FitTracks.Delete()



if float(h['DPangW'].Integral())!=0: Acc=float(h['DPangW'].Integral())/float(h['DPvesW'].Integral())
if float(h['DPangW_e'].Integral())!=0: Acc_e=float(h['DPangW_e'].Integral())/float(h['DPvesW_e'].Integral())
if float(h['DPangW_mu'].Integral())!=0: Acc_mu=float(h['DPangW_mu'].Integral())/float(h['DPvesW_mu'].Integral())
if float(h['DPangW_pi'].Integral())!=0: Acc_pi=float(h['DPangW_pi'].Integral())/float(h['DPvesW_pi'].Integral())
if float(h['DPangW_ka'].Integral())!=0: Acc_ka=float(h['DPangW_ka'].Integral())/float(h['DPvesW_ka'].Integral())
if float(h['DPangW_oth'].Integral())!=0: Acc_oth=float(h['DPangW_oth'].Integral())/float(h['DPvesW_oth'].Integral())

if float(h['DPangW'].Integral())!=0: RecEff=float(h['DPW'].Integral())/float(h['DPangW'].Integral())
if float(h['DPangW_e'].Integral())!=0: RecEff_e=float(h['DPW_e'].Integral())/float(h['DPangW_e'].Integral())
if float(h['DPangW_mu'].Integral())!=0: RecEff_mu=float(h['DPW_mu'].Integral())/float(h['DPangW_mu'].Integral())
if float(h['DPangW_pi'].Integral())!=0: RecEff_pi=float(h['DPW_pi'].Integral())/float(h['DPangW_pi'].Integral())
if float(h['DPangW_ka'].Integral())!=0: RecEff_ka=float(h['DPW_ka'].Integral())/float(h['DPangW_ka'].Integral())
if float(h['DPangW_oth'].Integral())!=0: RecEff_oth=float(h['DPW_oth'].Integral())/float(h['DPangW_oth'].Integral())

if float(h['DPang'].Integral())*2e20!=0: RecW=float(h['DPW'].Integral())/float(h['DPang'].Integral())*2e20
if float(h['DPang_e'].Integral())*2e20!=0: RecW_e=float(h['DPW_e'].Integral())/float(h['DPang_e'].Integral())*2e20
if float(h['DPang_mu'].Integral())*2e20!=0: RecW_mu=float(h['DPW_mu'].Integral())/float(h['DPang_mu'].Integral())*2e20
if float(h['DPang_pi'].Integral())*2e20!=0: RecW_pi=float(h['DPW_pi'].Integral())/float(h['DPang_pi'].Integral())*2e20
if float(h['DPang_ka'].Integral())*2e20!=0: RecW_ka=float(h['DPW_ka'].Integral())/float(h['DPang_ka'].Integral())*2e20
if float(h['DPang_oth'].Integral())*2e20!=0: RecW_oth=float(h['DPW_oth'].Integral())/float(h['DPang_oth'].Integral())*2e20


o1 = tmp1+"_rec.dat"
o2 = tmp1+"_e.dat"
o3 = tmp1+"_mu.dat"
o4 = tmp1+"_pi.dat"
o5 = tmp1+"_ka.dat"
o6 = tmp1+"_oth.dat"
o7 = tmp1+"_ves.dat"
o8 = tmp1+"_tab.dat"
o9 = tmp1+"_rves.dat"
a=open(o1,'w+')
b=open(o2,'w+')
c=open(o3,'w+')
d=open(o4,'w+')
e=open(o5,'w+')
f=open(o6,'w+')
g=open(o7,'w+')
kl=open(o8,'w+')
o9 = tmp1+"_rves.dat"
nl=open(o9,'w+')

a.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW, RecW, RecEff)) 
a.write('\n')
c.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_mu, RecW_mu, RecEff_mu))
c.write('\n')
d.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_pi, RecW_pi, RecEff_pi))
d.write('\n')
e.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_ka, RecW_ka, RecEff_ka)) 
e.write('\n')
f.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW_oth, RecW_oth, RecEff_oth))
f.write('\n')
g.write('%.4g %s %.8g %.8g %.8g %.8g %.8g %.8g' %(mass_mc, eps, Acc, Acc_e, Acc_mu, Acc_pi, Acc_ka, Acc_oth))
g.write('\n')
kl.write("%.4g %s %.8g %.8g %.8g %.8g " %(mass_mc, eps, nEvents, float(h['DPves'].Integral()), float(h['DPang'].Integral()), float(h['DP'].Integral())))
kl.write('\n')
nl.write('%.4g %s %.8g %.8g %.8g' %(mass_mc, eps, RecW, RecW, RecEff)) 
nl.write('\n')


a.close()
b.close()
c.close()
d.close()
e.close()
f.close()
g.close()
kl.close()
nl.close()

h['eff_P1'].Sumw2()
h['eff_Pt1'].Sumw2()
h['eff_Eta1'].Sumw2()
h['eff_P2'].Divide(h['eff_P1']*2e20, h['eff_P0'],1.,1.,"B")
h['eff_Pt2'].Divide(h['eff_Pt1']*2e20, h['eff_Pt0'],1.,1.,"B")
h['eff_Eta2'].Divide(h['eff_Eta1']*2e20,h['eff_Eta0'],1.,1.,"B")
h['eff_P2'].Draw("E")
h['eff_P2'].SetTitle("Momentum ;P(GeV/c);Eff.")
h['eff_Pt2'].Draw("E")
h['eff_Pt2'].SetTitle("Transverse Momentum ;P_{t}(GeV/c);Eff.")
h['eff_Eta2'].Draw("E")
h['eff_Eta2'].SetTitle("PseudoRapidity ;#eta;Eff.")	


ROOT.gStyle.SetOptStat(0)
hfile =tmp1+"_ana.root"
ROOT.gROOT.cd()
ut.writeHists(h,hfile)
