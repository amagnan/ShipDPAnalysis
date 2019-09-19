import ROOT as r

f0=open('kalanR.txt','r')

f1=open('kalanA.txt','r')

l_f0=f0.readlines() 

l_f1=f1.readlines()

def find(x_i,l_i):
    for i in l_i:
        i=i.replace('\n','')
        #print i
        if x_i==i:
            #print i
            return True
    return False

for i in l_f0:
    i=i.replace('\n','')
    t1=find(i,l_f1)
    #if not t1: print i[0], i[1], i[2]
    if not t1: print i
