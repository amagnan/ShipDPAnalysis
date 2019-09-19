import ROOT as r

#f0=open('Prod_r.txt','r')

#f0=open('mesuretim.txt','r')

#f1=open('ruretim.txt','r')

f0=open('s.txt','r')

f1=open('r.txt','r')


l_f0=f0.readlines() 

l_f1=f1.readlines()

def find(x_i,l_i):
    eosship =  r.gSystem.Getenv("EOSSHIP")  
    for i in l_i:
        i=i.replace('\n','')
        #print i
        #i=i+"_rec.root"
        if x_i==i:
            #print i
            #print x_i, i
            #f=r.TFile.Open(eosship+'/eos/experiment/ship/data/DarkPhoton/PBC-June-3/190826/reco/'+x_i)
            f=r.TFile.Open(eosship+'/eos/experiment/ship/data/DarkPhoton/PBC-June-3/190912/reco/'+x_i)
            try:
                sTree=f.cbmsim
                sTree.GetEntry()
                #print "sTree"
                fit=sTree.GetBranch("FitTracks")
                #print "fitTracks"
                fit.GetEntry()
                #f.Close()
                return True
            except:
                #f.Close()
                return False
    return False

for i in l_f0:
    i=i.replace('\n','')
    i=i.split(' ')
    if i[0]=='meson': x1 = i[0]+'_'+i[3]+'_mass'+i[1]+'_eps'+i[2]+'_rec.root'
    if not i[0]=='meson': x1 = i[0]+'_mass'+i[1]+'_eps'+i[2]+'_rec.root'
    #print x1
    t1=find(x1,l_f1)
    if not t1:
        if i[0]=='meson': print i[0], i[1], i[2], i[3]
        if not i[0]=='meson': print i[0], i[1], i[2]
