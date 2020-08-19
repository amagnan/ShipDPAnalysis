from array import array

for i in {'meson','pbrem','pbrem1','comb','comb1','qcd'}:
    for j in {'P','M'}:
        f=open("../data/200731/"+i+"_ErrorRate"+j+".csv","r")
        g=open("../data/200731/"+i+"_Error"+j+".csv","w")
        l=f.readlines()
        s=[]
        for x in l:
            x=x.split(",")
            g.write("%s,%s\n"%(x[0],str(float(x[1])**2.)))
        f.close()
        g.close()
