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

tmp=inputFile.replace('/eos/experiment/ship/data/DarkPhoton/PBC-June-3/sim/','')
tmp1=tmp.replace('.root','')
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
ut.bookHist(h,'DPvesW_oth','invariant Mass with Weights (GeV)',50,mass_mc-1.,mass_mc+1.)a

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

def isInFiducial(X,Y,Z):
	if Z > ShipGeo.TrackStation1.z : return False
	if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
	if dist2InnerWall(X,Y,Z)<5*u.cm: return False
	return True 
def angular(pdg,vtx,a):
	if

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
	LS = ROOT.gRandom.Uniform(flmin,flmax)
	p = dark.GetP()
	e = dark.GetEnergy()
	gam = e/ROOT.TMath.Sqrt(e*e-p*p)
	beta = p/e
	DP_instance = darkphoton.DarkPhoton(mass_mc,eps)
	ctau = DP_instance.cTau()
	return ROOT.TMath.Exp(-LS/(beta*gam*ctau))*( (ShipGeo.NuTauTarget.zdim*2.)/(beta*gam*ctau) )

def angle(X,Y,Z,Px,Py,Pz):
	x,y 	= X, Y
	m_x, m_y, m_z	= Px, Py, Pz
	a_x=np.arctan(m_x/m_z)
	a_y=np.arctan(m_y/m_z)
	return a_x,a_y

def invariant(S_x,S_y,e_0,e_1,p_0,p_1):
	a=r.TMath.Sqrt((S_x[0]-S_x[1])**2.+(S_y[0]-S_y[1])**2.)
	m_0=r.TMath.Sqrt(e_0**2-p_0**2)
	m_1=r.TMath.Sqrt(e_1**2-p_1**2)
	return r.TMath.Sqrt(m_0**2.+m_1**2.-2.*(e_0*e_1-p_0*p_1*m.cos(a)))

def myEventLoop(n,xsw):
	dp,	e,	mu,	tau,	pi,	ka,	oth	=0, 0, 0, 0, 0, 0, 0
	dp_a,	e_a,	mu_a,	tau_a,	pi_a, 	ka_a,	oth_a	=0, 0, 0, 0, 0, 0, 0
	dp_v, 	e_v, 	mu_v, 	tau_v	pi_v, 	ka_v,	oth_v	=0, 0, 0, 0, 0, 0, 0
	rc=sTree.GetEntry(n)
	had=[]
	mo=-99
	E,P,S_x,S_y=[],[],[],[]
	if pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,'meson',sTree.MCTrack[0].GetPdgCode())
	wg = sTree.MCTrack[1].GetWeight()
	if abs(wg)<0.00001: wg = wgnew(sTree.MCTrack[1])
	for mc,track in enumerate(sTree.MCTrack):
		if sTree.MCTrack[track.GetMotherId()].GetPdgCode()!=9900015 and sTree.MCTrack[track.GetMotherId()].GetPdgCode()!=4900023 and (track.GetPdgCode()==9900015 or track.GetPdgCode()==4900023):
			dp+=1
			mo=mc
			h['eff_P0w'].Fill(track.GetP(),wg*xsw)			
			h['eff_Pt0w'].Fill(track.GetPt(),wg*xsw)
			h['eff_Eta0w'].Fill(track.GetRapidity(),wg*xsw)

		if track.GetMotherId()==mo:
			h['eff_P0'].Fill(track.GetP())			
			h['eff_Pt0'].Fill(track.GetPt())
			h['eff_Eta0'].Fill(track.GetRapidity())
			if abs(track.GetPdgCode())==11 or abs(track.GetPdgCode())==13 or abs(track.GetPdgCode())==15:
				X=track.GetStartX()
				Y=track.GetStartY()
				Z=track.GetStartZ()
				Pz=track.GetPx()
				Py=track.GetPy()
				Pz=track.GetPz()
				P.append(track.GetP())
				E.append(track.GetEnergy())
				a=angle(X,Y,Z,Px,Py,Pz)
				S_x.append(a[0])
				S_y.append(a[1])
				if abs(track.GetPdgCode())==11:
					e+=1
					h['eff_P_e'].Fill(track.GetP())			
					h['eff_Pt_e'].Fill(track.GetPt())
					h['eff_Eta_e'].Fill(track.GetRapidity())
				if abs(track.GetPdgCode())==13:
					mu+=1
					h['eff_P_mu'].Fill(track.GetP())			
					h['eff_Pt_mu'].Fill(track.GetPt())
					h['eff_Eta_e'].Fill(track.GetRapidity())
				if abs(track.GetPdgCode())==15:
					tau+=1
					h['eff_P_tau'].Fill(track.GetP())			
					h['eff_Pt_tau'].Fill(track.GetPt())
					h['eff_Eta_tau'].Fill(track.GetRapidity())
			#else:#abs(1,2,3,4,5,6
			if abs(track.GetPdgCode())<7:
				had.append(mc)
			if abs(track.GetPdgCode())>6: print track.GetPdgCode()#for testing else state

		if len(had)==2:
			if track.GetMotherId()==had[0] or track.GetMotherId()==had[1]:
				X=track.GetStartX()
				Y=track.GetStartY()
				Z=track.GetStartZ()
				Pz=track.GetPx()
				Py=track.GetPy()
				Pz=track.GetPz()
				P.append(track.GetP())
				E.append(track.GetEnergy())
				a=angle(X,Y,Z,Px,Py,Pz)
				S_x.append(a[0])
				S_y.append(a[1])
				if abs(track.GetPdgCode())==211:
					pi+=1
					h['eff_P_pi'].Fill(track.GetP())			
					h['eff_Pt_pi'].Fill(track.GetPt())
					h['eff_Eta_pi'].Fill(track.GetRapidity())
				if abs(track.GetPdgCode())==321:
					ka+=1
					h['eff_P_ka'].Fill(track.GetP())			
					h['eff_Pt_ka'].Fill(track.GetPt())
					h['eff_Eta_ka'].Fill(track.GetRapidity())
				else:
					oth+=1
					h['eff_P_oth'].Fill(track.GetP())			
					h['eff_Pt_oth'].Fill(track.GetPt())
					h['eff_Eta_oth'].Fill(track.GetRapidity())
	if dp_v==2:
		h['eff_P2'].Fill(DP.P(),wg*xsw)
		h['eff_Pt2'].Fill(DP.Pt(),wg*xsw)
		h['eff_Eta2'].Fill(DP.Eta(),wg*xsw)
		h['DPves'].Fill(mass)
		h['DPvesW'].Fill(mass,wg*xsw)
		if e_v==2:  
			h['DPves_e'].Fill(mass)
			h['DPvesW_e'].Fill(mass,wg*xsw)
		if mu_v==2:              
			h['DPves_mu'].Fill(mass)
			h['DPvesW_mu'].Fill(mass,wg*xsw)
		if pi_v==2:              
			h['DPves_pi'].Fill(mass)
			h['DPvesW_pi'].Fill(mass,wg*xsw)
		if ka_v==2:
			h['DPves_ka'].Fill(mass)		
			h['DPvesW_ka'].Fill(mass,wg*xsw)

	if (e_v==1 or mu_v==1) and (pi_v==1 or ka_v==1):
		h['DPves_oth'].Fill(mass)
		h['DPvesW_oth'].Fill(mass,wg*xsw)

	if dp_a==2:
		h['eff_P1'].Fill(DP.P(),wg*xsw)
		h['eff_Pt1'].Fill(DP.Pt(),wg*xsw)
		h['eff_Eta1'].Fill(DP.Eta(),wg*xsw)
		h['DPang'].Fill(mass)
		h['DPangW'].Fill(mass,wg*xsw)
		if e_a==2:
			h['DPang_e'].Fill(mass)
			h['DPangW_e'].Fill(mass,wg*xsw)
		if mu_a==2:
			h['DPang_mu'].Fill(mass)
			h['DPangW_mu'].Fill(mass,wg*xsw)
		if pi_a==2:
			h['DPang_pi'].Fill(mass)
			h['DPangW_pi'].Fill(mass,wg*xsw)
		if ka_a==2:
			h['DPang_ka'].Fill(mass)
			h['DPangW_ka'].Fill(mass,wg*xsw)
	if (e_a==1 or mu_a==1) and (pi_a==1 or ka_a==1):
		oth_pair=1
		h['DPang_oth'].Fill(mass)
		h['DPangW_oth'].Fill(mass,wg*xsw)


nEvents =sTree.GetEntries()
print nEvents
for n in range(nEvents):
	if not pro=='meson': xsw = dputil.getDPprodRate(mass_mc,eps,pro,0)
	myEventLoop(n,xsw)
	sTree.FitTracks.Delete()



if float(h['DPvesW'].Integral())!=0: Acc=float(h['DPangW'].Integral())/float(h['DPvesW'].Integral())
if float(h['DPvesW_e'].Integral())!=0: Acc_e=float(h['DPangW_e'].Integral())/float(h['DPvesW_e'].Integral())
if float(h['DPvesW_mu'].Integral())!=0: Acc_mu=float(h['DPangW_mu'].Integral())/float(h['DPvesW_mu'].Integral())
if float(h['DPves_tau'].Integral())*2e20!=0: RecW_tau=float(h['DPangW_tau'].Integral())/float(h['DPves_tau'].Integral())*2e20
if float(h['DPvesW_pi'].Integral())!=0: Acc_pi=float(h['DPangW_pi'].Integral())/float(h['DPvesW_pi'].Integral())
if float(h['DPvesW_ka'].Integral())!=0: Acc_ka=float(h['DPangW_ka'].Integral())/float(h['DPvesW_ka'].Integral())
if float(h['DPvesW_oth'].Integral())!=0: Acc_oth=float(h['DPangW_oth'].Integral())/float(h['DPvesW_oth'].Integral())
if float(h['DPvesW_mix'].Integral())!=0: Acc_mix=float(h['DPangW_mix'].Integral())/float(h['DPvesW_mix'].Integral())

if float(h['DPves'].Integral())*2e20!=0: RecW=float(h['DPangW'].Integral())/float(h['DPves'].Integral())*2e20
if float(h['DPves_e'].Integral())*2e20!=0: RecW_e=float(h['DPangW_e'].Integral())/float(h['DPves_e'].Integral())*2e20
if float(h['DPves_mu'].Integral())*2e20!=0: RecW_mu=float(h['DPangW_mu'].Integral())/float(h['DPves_mu'].Integral())*2e20
if float(h['DPves_tau'].Integral())*2e20!=0: RecW_tau=float(h['DPangW_tau'].Integral())/float(h['DPves_tau'].Integral())*2e20
if float(h['DPves_pi'].Integral())*2e20!=0: RecW_pi=float(h['DPangW_pi'].Integral())/float(h['DPves_pi'].Integral())*2e20
if float(h['DPves_ka'].Integral())*2e20!=0: RecW_ka=float(h['DPangW_ka'].Integral())/float(h['DPves_ka'].Integral())*2e20
if float(h['DPves_oth'].Integral())*2e20!=0: RecW_oth=float(h['DPangW_oth'].Integral())/float(h['DPves_oth'].Integral())*2e20
if float(h['DPves_mix'].Integral())*2e20!=0: RecW_mix=float(h['DPangW_mix'].Integral())/float(h['DPves_mix'].Integral())*2e20


o1 = tmp1+"_rec.dat"
o2 = tmp1+"_e.dat"
o3 = tmp1+"_mu.dat"
o4 = tmp1+"_pi.dat"
o5 = tmp1+"_ka.dat"
o6 = tmp1+"_oth.dat"
o8 = tmp1+"_tab.dat"
a=open(o1,'w+')
b=open(o2,'w+')
c=open(o3,'w+')
d=open(o4,'w+')
e=open(o5,'w+')
f=open(o6,'w+')
kl=open(o8,'w+')



a.write('%.4g %s %.8g %.8g %.8g' %(mass, eps, RecW, RecW, Acc)) 
a.write('\n')
b.write('%.4g %s %.8g %.8g %.8g' %(mass, eps, RecW_e, RecW_e, Acc_e))
b.write('\n')
c.write('%.4g %s %.8g %.8g %.8g' %(mass, eps, RecW_mu, RecW_mu, Acc_mu))
c.write('\n')
d.write('%.4g %s %.8g %.8g %.8g' %(mass, eps, RecW_pi, RecW_pi, Acc_pi))
d.write('\n')
e.write('%.4g %s %.8g %.8g %.8g' %(mass, eps, RecW_ka, RecW_ka, Acc_ka)) 
e.write('\n')
f.write('%.4g %s %.8g %.8g %.8g' %(mass, eps, RecW_oth, RecW_oth, Acc_oth))
f.write('\n')
kl.write("%.4g %s %.8g %.8g %.8g %.8g " %(mass, eps, nEvents, float(h['DP'].Integral()), float(h['DPves'].Integral()), float(h['DPang'].Integral())))
kl.write('\n')


a.close()
b.close()
c.close()
d.close()
e.close()
f.close()
kl.close()


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
