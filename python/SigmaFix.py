date='200826'
#date='alp'
import math
#prods=['comb','comb1','meson','meson_pi0','meson_omega','meson_eta','meson_eta1','pbrem','pbrem1','qcd']
#prods=['meson','meson_pi0','meson_omega','meson_eta','meson_eta1','meson_eta11']
prods=['qcd']
#prods=['ALPACA']
for prod in prods:
    #exec('f1 = open("../data/"+date+"/%s"+"_FixedRate1.dat","w")'%(prod))
    f1 = open("../data/"+date+"/"+prod+"_Rate1.dat","w")
    exec('f0_%s = open("../data/"+date+"/%s"+"_rate1.dat","r")'%(prod,prod))
    exec('l0 = f0_%s.readlines() '%(prod))
    for i in l0:
        i=i.replace('\n','')
        i=i.split(' ')
        #if not ("meson" in prod): i[2]=str(float(i[2])*3.7817749617068297)
        #else: i[2]=str(float(i[2])/0.6)
        i[2]=str(float(i[2])*3.7817749617068297)
        #i[2]=str(float(i[2])*0.6)
        #print type(i[2]),i[2]
        f1.write("%s %s %s"%(i[0],i[1],i[2]))
        f1.write("\n")
    f1.close()
