import pandas as pd
try:
    opts, args = getopt.getopt(sys.argv[1:], "d:", ["date="]) 
except getopt.GetoptError:
    sys.exit()

for o,a in opts:
    if o in ('-d','--date',): date = a

prods = ['pbrem','pbrem1','meson','qcd','comb','comb1']
for prod in prods:
    inp="../data/"+date+"/"+prod+"_Rate1.csv"
    data = pd.read_csv(path,header=None)
    data.head()
    #new=data.sort_values(by=[1,0],ascending=[False,True])
    new=data.sort_values(by=1,ascending=False)
    new.to_csv(inp,index=False) 

