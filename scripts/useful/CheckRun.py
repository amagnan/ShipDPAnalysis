import ROOT as r

n=raw_input()

f0=open('mesuretim.txt','r')

f1=open('muretim'+n+'.txt','r')

l_f0=f0.readlines() 

l_f1=f1.readlines()

def find(x_i,l_i):
    eosship =  r.gSystem.Getenv("EOSSHIP")
    for i in l_i:
        i=i.replace('\n','')
        #print i
        if x_i==i:
            f=r.TFile.Open(eosship+'/eos/experiment/ship/data/DarkPhoton/PBC-June-3/190826/runs/'+x_i)
            try:
                sTree=f.cbmsim
                #print "sTree"
                #fit=sTree.GetBranch("FitTracks")
                #print "fitTracks"
                sTree.GetEntry()
                #f.Close()
                return True
            except:
                #f.Close()
                return False
    return False

for i in l_f0:
    i=i.replace('\n','')
    i=i.split(' ')
    x1=i[0]+'_'+i[3]+'_mass'+i[1]+'_eps'+i[2]+'_'+n+'.root'
    t1=find(x1,l_f1)
    if not t1: print i[0], i[1], i[2], n, i[3]
