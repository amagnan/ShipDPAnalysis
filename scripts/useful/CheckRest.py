import ROOT as r

f0=open('yeniuretim.txt','r')

f1=open('sonuretim.txt','r')

l_f0=f0.readlines() 

l_f1=f1.readlines()

eosship =  r.gSystem.Getenv("EOSSHIP")

def findN(x_i,l_i):
    eosship =  r.gSystem.Getenv("EOSSHIP")
    for i in l_i:
        i=i.replace('\n','') 
        #print i
        if x_i==i:
            #print i
            f=r.TFile.Open(eosship+'/eos/experiment/ship/data/DarkPhoton/PBC-June-3/190420/sim/'+x_i)
            if f:
                k=f.GetListOfKeys()
                if k.Contains("cbmsim"):
                    sTree=f.cbmsim
                    nEvents=sTree.GetEntries()
                    #print nEvents
                    if abs( nEvents-5000.)/5000.<0.11: 
                        #sTree.Close()
                        f.Close()
                        return True
                    else:
                        print nEvents
                        #sTree.Close()
                        f.Close()
                        return False
                else:
                    #sTree.Close()
                    f.Close()
                    return False
            else: 
                #f.close()
                return False
            #return True
    return False

for i in l_f0:
    i=i.replace('\n','')
    i=i.split(' ')
    x1=i[0]+'_mass'+i[1]+'_eps'+i[2]+'.root'
    a1=findN(x1,l_f1)
    if not a1: print i[0], i[1], i[2]
