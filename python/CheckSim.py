import ROOT as r

f0=open('New.txt','r')

f1=open('Prod_1.txt','r')
f2=open('Prod_2.txt','r')
f3=open('Prod_3.txt','r')
f4=open('Prod_4.txt','r')
f5=open('Prod_5.txt','r')

l_f0=f0.readlines() 

l_f1=f1.readlines()
l_f2=f2.readlines()
l_f3=f3.readlines()
l_f4=f4.readlines()
l_f5=f5.readlines()

eosship =  r.gSystem.Getenv("EOSSHIP")

def findN(x_i,l_i):
    eosship =  r.gSystem.Getenv("EOSSHIP")
    for i in l_i:
        i=i.replace('\n','') 
        #print i
        if x_i==i:
            #print i
            f=r.TFile.Open(eosship+'/eos/experiment/ship/data/DarkPhoton/PBC-June-3/190305/runs/'+x_i)
            if f:
                k=f.GetListOfKeys()
                if k.Contains("cbmsim"):
                    sTree=f.cbmsim
                    nEvents=sTree.GetEntries()
                    if abs( nEvents-1000.)/1000.<0.11: 
                        #sTree.Close()
                        f.Close()
                        return True
                    else:
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
    return False

for i in l_f0:
    i=i.replace('\n','')
    i=i.split(' ')
    x1=i[0]+'_mass'+i[1]+'_eps'+i[2]+'_1.root'
    x2=i[0]+'_mass'+i[1]+'_eps'+i[2]+'_2.root'
    x3=i[0]+'_mass'+i[1]+'_eps'+i[2]+'_3.root'
    x4=i[0]+'_mass'+i[1]+'_eps'+i[2]+'_4.root'
    x5=i[0]+'_mass'+i[1]+'_eps'+i[2]+'_5.root'
    #if float(i[1])<0.30: continue
    a1=findN(x1,l_f1)
    if not a1: print i[0], i[1], i[2], "1" 
    a2=findN(x2,l_f2)
    if not a2: print i[0], i[1], i[2], "2" 
    a3=findN(x3,l_f3)
    if not a3: print i[0], i[1], i[2], "3" 
    a4=findN(x4,l_f4)
    if not a4: print i[0], i[1], i[2], "4" 
    a5=findN(x5,l_f5)
    if not a5: print i[0], i[1], i[2], "5"
