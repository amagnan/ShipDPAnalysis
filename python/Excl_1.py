import pandas as pd
import numpy as np
import os,sys,getopt
from array import array

leptophilic=0 

try:
    opts, args = getopt.getopt(sys.argv[1:], "l:")

except getopt.GetoptError:
    print 'decay to all SM particles'

for o,a in opts:
    if o in ('-l',): leptophilic = a 

pathR = "../Exclusion/"

if leptophilic:
    data1 = pd.read_csv(pathR+'qcd_Rate2.csv', header=None)
    data2 = pd.read_csv(pathR+'meson_Rate2.csv', header=None)
    data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
    dataC = pd.read_csv(pathR+'comb_Rate2.csv', header=None)
else:
    data1 = pd.read_csv(pathR+'qcd_Rate1.csv', header=None)
    data2 = pd.read_csv(pathR+'meson_Rate1.csv', header=None)
    data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
    dataC = pd.read_csv(pathR+'comb_Rate1.csv', header=None)
    #dataC = pd.read_csv(pathW+'comb.csv', header=None)
##### Constraints
dataC1 = pd.read_csv(pathR+'Dataset_DP_excluded_u.csv', header=None)
dataC2 = pd.read_csv(pathR+'Dataset_DP_excluded_l.csv', header=None)

y2 = np.array([1e-2]*len(dataC1[1]))
#y3=np.array([1.0e-4]*len(data3[1]))
#y4=np.array([200.00]*len(data3[1]))

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as pl

pl.yscale('log')
pl.xscale('log')

pl.plot(data1[0], data1[1], 'b', linewidth=2.0, label='SHiP (QCD)', labelsize=20)
pl.plot(data3[0], data3[1], 'g', linewidth=2.0,label='SHiP (Pbrem)', labelsize=20)
pl.plot(data2[0], data2[1], 'r', linewidth=2.0,label='SHiP (Mesons)', labelsize=20)
pl.plot(dataC[0], dataC[1], 'black', linewidth=1.8, label='SHiP (Combined)', labelsize=20)

#pl.fill_between(data1[0], data1[1], y2=0, color='red', hatch='.', label='SHiP (QCD)')
#pl.fill_between(data3[0], data3[1], y2=0, color='blue', label="SHiP (Brem)") 
#pl.fill_between(data2[0], data2[1], y2=0, color='green', label="SHiP (Mesons)")

#pl.fill_between(data3[0], data3[1], where=data3[0] > y3  , color='green', label="SHiP (Brem)")

pl.fill_between(dataC1[0], dataC1[1], y2=y2, color="grey", label="Excluded Region", labelsize=20)
pl.fill_between(dataC2[0], dataC2[1], y2=0, color="grey")
 
pl.xlabel("M [MeV]",fontsize=20)
pl.ylim(1e-10, 1e-2)
pl.xlim(13, 1e4)
pl.ylabel(r"$\epsilon$", fontsize=20)
pl.legend(prop={'size': 10})
pl.savefig(pathR+"excl_"+str(leptophilic)+".pdf")
pl.show()
