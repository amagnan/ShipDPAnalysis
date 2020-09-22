date='200826'
import math
#date='ALPACA/200820'
#prods=['meson','pbrem','pbrem1','qcd','comb','comb1']
prods=['comb','comb1','qcd']
modes=['Rate1','ErrorRateP','ErrorRateM']
#prods=['ALPACA']
#modes =['rate']
for mode in modes:
    for prod in prods:
        f0 = open("../data/"+date+"/"+prod+"_"+mode+".csv","r")
        f1 = open("../data/"+date+"/"+prod+"_"+mode+"_Ordered.csv","w")
        l0 = f0.readlines()
        i=0
        n=len(l0)
        k=[['none']*2]*n
        for l in l0:
            l=l.replace('\n','')
            l=l.split(', ')
            if i%2==0:
                k[i/2]=[str(l[0]),str(l[1])]
            if i%2!=0:
                k[n-(i+1)/2]=[str(l[0]),str(l[1])]
            i+=1
        for m in k:
            f1.write("%s, %s"%(str(m[0]),str(float(m[1])**2.)))
            #f1.write("%s, %s"%(str(m[0]),str(float(m[1]))))
            f1.write("\n")
        f1.close()
