import ROOT as r

#n=raw_input()

f0=open('uretim.txt','r')
f1=open('mesuretim.txt','r')

l_f0=f0.readlines() 
l_f1=f1.readlines()

def find(x_i,l_i):
    for i in l_i:
        i=i.replace('\n','')
        i=i.split(' ')
        i=i[0]+'_'+i[3]+'_mass'+i[1]+'_eps'+i[2]
        #print i
        if x_i==i:
            return True
    return False

for i in l_f0:
    i=i.replace('\n','')
    i=i.split(' ')
    x1=i[0]+'_'+i[3]+'_mass'+i[1]+'_eps'+i[2]
    t1=find(x1,l_f1)
    if not t1: print x1
