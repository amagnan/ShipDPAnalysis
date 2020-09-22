import pandas as pd
import os,sys,getopt
try:
    opts, args = getopt.getopt(sys.argv[1:], "d:", ["date="]) 
except getopt.GetoptError:
    sys.exit()

for o,a in opts:
    if o in ('-d','--date',): date = a

#prods = ['pbrem','pbrem1','meson','qcd','comb','comb1']
prods = ['qcd']
#for i in {"Rate1","ErrorRateM","ErrorRateP"}:
for i in {"Rate"}:
    for prod in prods:
        inp="../data/"+date+"/"+prod+"_"+i+".dat"
        data = pd.read_csv(inp,header=None,sep="\s+")
        #data.head()
        #new=data.sort_values(by=[1,0],ascending=[False,True])
        new=data.sort_values(by=0,ascending=True)
        new.to_csv(inp,index=False,header=None) 
