import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
import proton_bremsstrahlung
import dpProductionRates as dputil
import darkphoton
import sys

import shipRoot_conf
shipRoot_conf.configure()


debug = False
PDG = ROOT.TDatabasePDG.Instance()
inputFile  = None
geoFile    = None
nEvents    = 9999999
protonFlux = 2e20
mass = 0.5
epsilon = 3e-8
prodMode = "pbrem"

flmin=47.505*u.m
flmax=98.265*u.m
ROOT.gRandom.SetSeed(0)

zmin=-2477.5#cm
#(ship_geo.Chamber1.z - ship_geo.chambers.Tub1length) - ship_geo.target.z0
zmax=2598.5#cm
#ship_geo.TrackStation1.z - ship_geo.target.z0
xmax=2.5*u.m
ymax=5*u.m


try:
    opts, args = getopt.getopt(sys.argv[1:], "n:f:g:m:e:p", ["nEvents=","geoFile=","epsilon=","prodmode="])
except (getopt.GetoptError or len(args<2)):
    # print help information and exit:
    print 'Usage: %s -f <filename>'%sys.argv[0]
    sys.exit()
for o, a in opts:
    if o in ("-f",):
        inputFile = a
    if o in ("-m",):
        mass = float(a)
    if o in ("-e","--epsilon"):
        epsilon = float(a)
    if o in ("-g", "--geoFile",):
        geoFile = a
    if o in ("-n", "--nEvents=",):
        nEvents = int(a)
    if o in ("-p","--prodmode",):
        prodMode = a

eosship = ROOT.gSystem.Getenv("EOSSHIP")
if inputFile[0:4] == "/eos":
  eospath = eosship+inputFile
  f = ROOT.TFile.Open(eospath)
  sTree = f.cbmsim
else:
  f = ROOT.TFile(inputFile)
  sTree = f.cbmsim


h = {}
ut.bookHist(h,'DPgenM','generated Mass (GeV)',120,0.,12.)
ut.bookHist(h,'DPgenMw','generated Mass with weights (GeV)',120,0.,12.)
ut.bookHist(h,'DPvtxMw','generated Mass with weights (GeV)',120,0.,12.)
ut.bookHist(h,'DPselMw','selected Mass with weights (GeV)',120,0.,12.)
ut.bookHist(h,'DPselMwnew','selected Mass with new weights (GeV)',120,0.,12.)
ut.bookHist(h,'rate','Total production rate (events/2e20 p.o.t.)',1000,0,5000)
ut.bookHist(h,'rateNew','Total production rate (events/2e20 p.o.t.)',1000,0,5000)
ut.bookHist(h,'DPdecay','Decay pdg',1200,-600,600)
ut.bookHist(h,'eEta','Decay e #eta',100,0,10)
ut.bookHist(h,'muEta','Decay #mu #eta',100,0,10)
ut.bookHist(h,'tauEta','Decay #tau #eta',100,0,10)
ut.bookHist(h,'hadEta','Decay had #eta',100,0,10)
ut.bookHist(h,'ePT','Decay e p_{T} (GeV)',1000,0,50)
ut.bookHist(h,'muPT','Decay #mu p_{T} (GeV)',1000,0,50)
ut.bookHist(h,'tauPT','Decay #tau p_{T} (GeV)',1000,0,50)
ut.bookHist(h,'hadPT','Decay had p_{T} (GeV)',1000,0,50)
ut.bookHist(h,'ePhi','Decay e #phi',100,-3.1416,3.1416)
ut.bookHist(h,'muPhi','Decay #mu #phi',100,-3.1416,3.1416)
ut.bookHist(h,'tauPhi','Decay #tau #phi',100,-3.1416,3.1416)
ut.bookHist(h,'hadPhi','Decay had #phi',100,-3.1416,3.1416)

ut.bookHist(h,'eeM','Decay ee M (GeV)',1000,0,10)
ut.bookHist(h,'mumuM','Decay #mu#mu M (GeV)',1000,0,10)
ut.bookHist(h,'tautauM','Decay #tau#tau M (GeV)',1000,0,10)
ut.bookHist(h,'hadM','Decay had M (GeV)',1000,0,10)


def makePlots():
    ut.bookCanvas(h,key='DPanalysis',title='Mass',nx=800,ny=600,cx=2,cy=1)
    cv = h['DPanalysis'].cd(1)
    h['DPgenM'].Draw()
    cv = h['DPanalysis'].cd(2)
    h['DPselMw'].Draw()
    h['DPanalysis'].Print("./DPmass.pdf")
    h['DPanalysis'].Print("./DPmass.png")
    print 'finished making plots'

def checkMass(track):
    mymass = ROOT.TMath.Sqrt(sTree.MCTrack[track].GetEnergy()**2-sTree.MCTrack[track].GetP()**2)
    if (abs(mymass-mass)>0.001):
        print "Warning! Mass in file %3.3f GeV not same as mass given in argument: %3.3f GeV"%(mymass,mass)
        sys.exit()
    return mymass

def findDP():
    select=-1
    for dpidx in range(0,sTree.MCTrack.GetEntries()):
        dpid = sTree.MCTrack[dpidx].GetPdgCode()
        if dputil.isDP(dpid): 
            select=dpidx
            mymass = checkMass(select)
            #print 'selected %d with pdgid %d, mass %3.3f GeV'%(select,dpid,mymass)
    return select

def findMumPDGID():
    select=0
    for dpidx in range(0,sTree.MCTrack.GetEntries()):
        dpid = sTree.MCTrack[dpidx].GetPdgCode()
        #print "%d %d %d"%(dpidx,dpid,mum)
        if dputil.isDP(dpid):
            mum = sTree.MCTrack[dpidx].GetMotherId()
            #print 'selected mum %d'%mum
            select = sTree.MCTrack[mum].GetPdgCode()
            mymass = checkMass(dpidx)
            #print 'selected mum pdgid %d for pdgid %d with mass %3.3f'%(select,dpid,mymass)
            break
    return select

def findDecayProd(dpidx):
    daughters=[]
    #print "DP index: %d"%dpidx
    #for leptons: retrieve the indices of the last particles in the list with the right flavour.
    # For hadrons: direct decay products of the DP.
    for idx in range(0,sTree.MCTrack.GetEntries()):
        mum = sTree.MCTrack[idx].GetMotherId()
        #print "mum: %d"%mum
        if (mum==dpidx):
            #print "%d %d %d %3.3f %3.3f"%(idx,mum,sTree.MCTrack[idx].GetPdgCode(),sTree.MCTrack[idx].GetPt(),sTree.MCTrack[idx].GetRapidity())
            pdg=sTree.MCTrack[idx].GetPdgCode()
            if (abs(pdg)>10 and abs(pdg) < 16):
                daughters.append(idx)
            else:
                if (pdg!=22):
                    daughters.append(idx)

#rescue second lepton
    if (len(daughters)==1):
        pdg=sTree.MCTrack[daughters[0]].GetPdgCode()
        if (abs(pdg)>10 and abs(pdg) < 16):
            for idx in range(0,sTree.MCTrack.GetEntries()):
                if (sTree.MCTrack[idx].GetPdgCode()==-pdg):
                    daughters.append(idx)
                    break
    
    return daughters

# start event loop
def myEventLoop(n):
    if (n%1000==0): 
        print "--Processing event %d"%n
    rc = sTree.GetEntry(n)
    nTrks = sTree.MCTrack.GetEntries()
    dpindex = findDP()
    counts=[0,0,0,0,0,0,0,0]
    if (dpindex < 0):
        print "ERROR! Did not find dark photon in event %d..."%n
        counts[0] += 1
        return counts


    if (prodMode=="meson"):
        xsw = dputil.getDPprodRate(mass,epsilon,prodMode,findMumPDGID())
    else:
        xsw = dputil.getDPprodRate(mass,epsilon,prodMode,0)


    if (xsw == 0):
        print "ERROR! Production rate is 0 in event %d..."%n
        counts[1] += 1

    mymass = checkMass(dpindex)

    wg = sTree.MCTrack[dpindex].GetWeight()
    LS = ROOT.gRandom.Uniform(flmin,flmax)
    p = sTree.MCTrack[dpindex].GetP()
    e = sTree.MCTrack[dpindex].GetEnergy()
    gam  = e/ROOT.TMath.Sqrt(e*e-p*p)
    beta = p/e 
    DP_instance = darkphoton.DarkPhoton(mass,epsilon)
    ctau = DP_instance.cTau()
    wnew = ROOT.TMath.Exp(-LS/(beta*gam*ctau))*( (flmax-flmin)/(beta*gam*ctau) )

    #print "A' decay vertex rescaling weight: %.8g %.8g %.8g %.8g %.8g"%(beta,gam,ctau,wg,wnew)
    #print "A' info: %3.3f %d %d %3.3f %3.3f" %(mass,dpindex,sTree.MCTrack[dpindex].GetPdgCode(),sTree.MCTrack[dpindex].GetPt(),sTree.MCTrack[dpindex].GetRapidity())
    
    #select only decays to mu
    daughters=findDecayProd(dpindex)
    for p in range(0,len(daughters)):
        h['DPdecay'].Fill(sTree.MCTrack[daughters[p]].GetPdgCode())
    if (len(daughters)!=2):
        counts[6] += 1
        print "Event %d, Daughters:\t"%n
        for p in range(0,len(daughters)):
            print "pdg %d \t"%(sTree.MCTrack[daughters[p]].GetPdgCode())
        #for idx in range(0,sTree.MCTrack.GetEntries()):
            #print "%d %d %d %3.3f %3.3f"%(idx,sTree.MCTrack[idx].GetMotherId(),sTree.MCTrack[idx].GetPdgCode(),sTree.MCTrack[idx].GetPt(),sTree.MCTrack[idx].GetRapidity())
    else:
        h['DPgenM'].Fill(mymass)
        h['DPgenMw'].Fill(mymass,wnew*xsw)

        #filter with vertex in decay volume:
        decVtxx=sTree.MCTrack[daughters[0]].GetStartX()
        decVtxy=sTree.MCTrack[daughters[0]].GetStartY()
        decVtxz=sTree.MCTrack[daughters[0]].GetStartZ()
        if (decVtxz < zmin or decVtxz > zmax or abs(decVtxx)>xmax or abs(decVtxy)>ymax):
            counts[7] += 1
            return counts

        h['DPvtxMw'].Fill(mymass,wnew*xsw)

        pt1=sTree.MCTrack[daughters[0]].GetPt()
        pt2=sTree.MCTrack[daughters[1]].GetPt()
        eta1=ROOT.TMath.ASinH(sTree.MCTrack[daughters[0]].GetPz()/sTree.MCTrack[daughters[0]].GetPt())
        eta2=ROOT.TMath.ASinH(sTree.MCTrack[daughters[1]].GetPz()/sTree.MCTrack[daughters[1]].GetPt())
        px1 = sTree.MCTrack[daughters[0]].GetPx()
        py1 = sTree.MCTrack[daughters[0]].GetPy()
        phi1 = ROOT.TMath.ATan(py1/px1)
        if (px1<0): phi1 += abs(py1)/py1*3.1416
        px2 = sTree.MCTrack[daughters[1]].GetPx()
        py2 = sTree.MCTrack[daughters[1]].GetPy()
        phi2 = ROOT.TMath.ATan(py2/px2)
        if (px2<0): phi2 += abs(py2)/py2*3.1416


        lv1=ROOT.TLorentzVector(sTree.MCTrack[daughters[0]].GetPx(),sTree.MCTrack[daughters[0]].GetPy(),sTree.MCTrack[daughters[0]].GetPz(),sTree.MCTrack[daughters[0]].GetEnergy())
        lv2=ROOT.TLorentzVector(sTree.MCTrack[daughters[1]].GetPx(),sTree.MCTrack[daughters[1]].GetPy(),sTree.MCTrack[daughters[1]].GetPz(),sTree.MCTrack[daughters[1]].GetEnergy())
        dauMass = (lv1+lv2).M()
        #dauMass = 2*pt1*pt2*(ROOT.TMath.CosH(eta1-eta2)-ROOT.TMath.Cos(phi1-phi2))
        if (abs(sTree.MCTrack[daughters[0]].GetPdgCode())==13 and
            abs(sTree.MCTrack[daughters[1]].GetPdgCode())==13):
        #print "Select decay to %s %s"%(sTree.MCTrack[daughters[0]].GetPdgCode(),sTree.MCTrack[daughters[1]].GetPdgCode())
            h['DPselMw'].Fill(mymass,wg*xsw)
            h['DPselMwnew'].Fill(mymass,wnew*xsw)
            counts[3] += 1
            h['muEta'].Fill(eta1,wnew)
            h['muPT'].Fill(pt1,wnew)
            h['muPhi'].Fill(phi1,wnew)
            h['muEta'].Fill(eta2,wnew)
            h['muPT'].Fill(pt2,wnew)
            h['muPhi'].Fill(phi2,wnew)
            h['mumuM'].Fill(dauMass,wnew*xsw)
        else:
        #select only decays to e
            if (abs(sTree.MCTrack[daughters[0]].GetPdgCode())==11 and
                abs(sTree.MCTrack[daughters[1]].GetPdgCode())==11):
            #print "Reject decay to %s %s"%(sTree.MCTrack[daughters[0]].GetPdgCode(),sTree.MCTrack[daughters[1]].GetPdgCode())
                counts[2] += 1
                h['eEta'].Fill(eta1,wnew)
                h['ePT'].Fill(pt1,wnew)
                h['ePhi'].Fill(phi1,wnew)
                h['eEta'].Fill(eta2,wnew)
                h['ePT'].Fill(pt2,wnew)
                h['ePhi'].Fill(phi2,wnew)
                h['eeM'].Fill(dauMass,wnew*xsw)
        #select only decays to tau
            elif (abs(sTree.MCTrack[daughters[0]].GetPdgCode())==15 and
                  abs(sTree.MCTrack[daughters[1]].GetPdgCode())==15):
            #print "Reject decay to %s %s"%(sTree.MCTrack[daughters[0]].GetPdgCode(),sTree.MCTrack[daughters[1]].GetPdgCode())
                counts[4] += 1
                h['tauEta'].Fill(eta1,wnew)
                h['tauPT'].Fill(pt1,wnew)
                h['tauPhi'].Fill(phi1,wnew)
                h['tauEta'].Fill(eta2,wnew)
                h['tauPT'].Fill(pt2,wnew)
                h['tauPhi'].Fill(phi2,wnew)
                h['tautauM'].Fill(dauMass,wnew*xsw)
        #select only decays to hadrons
            else:
                #print "Reject decay to %s %s in event %d"%(sTree.MCTrack[daughters[0]].GetPdgCode(),sTree.MCTrack[daughters[1]].GetPdgCode(),n)
                counts[5] += 1
                h['hadEta'].Fill(eta1,wnew)
                h['hadPT'].Fill(pt1,wnew)
                h['hadPhi'].Fill(phi1,wnew)
                h['hadEta'].Fill(eta2,wnew)
                h['hadPT'].Fill(pt2,wnew)
                h['hadPhi'].Fill(phi2,wnew)
                h['hadM'].Fill(dauMass,wnew*xsw)


    return counts
    

#main action....

if (nEvents > 0):
    nEvents = min(sTree.GetEntries(),nEvents)
else:
    nEvents = sTree.GetEntries()



#debug info...
#if prodMode=="meson": 
#    dputil.getDPprodRate(prodMode,111,True)
#    dputil.getDPprodRate(prodMode,221,True)
#    dputil.getDPprodRate(prodMode,223,True)
#    dputil.getDPprodRate(prodMode,331,True)
#else: dputil.getDPprodRate(prodMode,0,True)

#for tries in range(0,0):
    #print " - try #%d"%tries

counts = [0,0,0,0,0,0,0,0]
print "-- Processing %d entries"%nEvents
sys.stdout.flush()
for n in range(0,nEvents): 
#for n in range(0,100): 
    tmpcounts = myEventLoop(n);
    counts[0] += tmpcounts[0]
    counts[1] += tmpcounts[1]
    counts[2] += tmpcounts[2]
    counts[3] += tmpcounts[3]
    counts[4] += tmpcounts[4]
    counts[5] += tmpcounts[5]
    counts[6] += tmpcounts[6]
    counts[7] += tmpcounts[7]
    sys.stdout.flush()
   
print " -- Info: events with 0 DP or rate = 0 or decay to 2 electrons, 2 muons, 2 taus, 2 hadrons, or others (not 2 daughters), vtx not in vessel: %s. Sum: %d"%(counts,sum(counts))
    
if h['DPgenM'].Integral() != 0:
    print "Ngen = %d"%(h['DPgenM'].Integral())
    if h['DPgenMw'].Integral() != 0:
        vtxAcc = h['DPvtxMw'].Integral()/h['DPgenMw'].Integral()
    else:
        vtxAcc=1
    
    decayW =  h['DPselMw'].Integral()/h['DPgenM'].Integral()
    decayWnew =  h['DPselMwnew'].Integral()/h['DPgenM'].Integral()
    print "Decay vertex in vessel acceptance: %.8g"%vtxAcc
    print "Decay vertex rescaling average weight:\t %.8g \t %.8g"%(decayW,decayWnew)
    print "Proton flux: %.8g"%protonFlux
    print "Total A' production rate: %s \t %.2g \t %.8g \t %.8g \t %.8g"%(prodMode,mass,epsilon,protonFlux*decayW,protonFlux*decayWnew)
    h['rate'].Fill(protonFlux*decayW)
    h['rateNew'].Fill(protonFlux*decayWnew)
else: 
    print " --- Problem, histograms empty, cannot get production rate!"
    print " --- h['DPselMw'].Integral() = %.8g, h['DPselMwnew'].Integral() = %.8g, h['DPgenM'].Integral() = %.8g"%(h['DPselMw'].Integral(),h['DPselMwnew'].Integral(),h['DPgenM'].Integral())


#makePlots()

tmp = inputFile.split('/')
#if (len(tmp)==1): hfile = "testAna_%s.root"%prodMode
#else : hfile = tmp[len(tmp)-1] 
hfile = "histAna_%s_%.2gGeV_%.8g.root"%(prodMode,mass,epsilon)
ROOT.gROOT.cd()
ut.writeHists(h,hfile)
