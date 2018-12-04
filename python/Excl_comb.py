import os,sys,getopt
 
from array import array
import pandas as pd
import numpy as np

leptophilic=0 

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:l:")

except getopt.GetoptError:
    print 'decay to all SM particles'

for o,a in opts:
    if o in ('-l',): leptophilic = a
    if o in ('-d',): date = a

pathW = "../data/"+date+"/"
pathR = "../Exclusion/"
##### Signal Datasets
 
if leptophilic:
    data1 = pd.read_csv(pathW+'combined_Ana_rate2.csv', header=None)
else:
    data1 = pd.read_csv(pathW+'combined_Ana_rate1.csv', header=None)
 
#### Constraints
dataC1 = pd.read_csv(pathR+'Dataset_DP_excluded_u.csv', header=None)
dataC2 = pd.read_csv(pathR+'Dataset_DP_excluded_l.csv', header=None)


y2 = np.array([1e-2]*len(dataC1[1]))
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as pl


pl.yscale('log')
pl.xscale('log')

pl.fill_between(data1[0], data1[1], y2=0, color='red', label="SHiP")
pl.fill_between(dataC1[0], dataC1[1], y2=y2, color="grey", label="Excluded Region")
pl.fill_between(dataC2[0], dataC2[1], y2=0, color="grey")
pl.xlabel("M [MeV]",fontsize=18)
pl.ylim(1e-8, 1e-2)
pl.xlim(13, 1e4)
pl.ylabel(r"$\epsilon$", fontsize=18)
pl.legend(prop={'size': 10})
pl.savefig("Excl_comb.pdf")
pl.show()
