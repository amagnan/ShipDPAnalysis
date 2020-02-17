date='200213'
import math
prod='pbrem'
mothers=['e','mu','tau','charged']
for mother in mothers:
    exec('f1_%s = open("../data/"+date+"/"+prod+"_%s"+"_Rate1.dat","r")'%(mother,mother))
    exec('l1_%s = f1_%s.readlines() '%(mother,mother))

f0 = open(prod+"_mass.txt","r")#dosya lazim
l0 = f0.readlines()

f1 = open("../data/"+date+"/"+prod+"_allCharged.dat","w")

def find(l0,l_i):
    for i in l_i:
        i=i.replace('\n','')
        i = i.split(' ')
        if float(l0[0])==float(i[0]) and float(l0[1])==float(i[1]): return float(i[2]),float(i[3])
    return 0.


for i in l0:
    s1,s2 = 0.,0.
    i=i.replace('\n','')
    i=i.split(' ')
    for mother in mothers:
        exec('s1 += find(i,l1_%s)'%(mother))
    if s1:
        f1.write('%s %s %s'%(i[0],i[1],s1[0],s1[1]))
        f1.write('\n')
