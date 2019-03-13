import pandas as pd
import numpy as np
import os,sys,getopt
from array import array

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

if leptophilic:
    data1 = pd.read_csv(pathW+'qcd_Ana_rate2.csv', header=None)
    data2 = pd.read_csv(pathW+'meson_Ana_rate2.csv', header=None)
    data3 = pd.read_csv(pathW+'pbrem_Ana_rate2.csv', header=None)
    dataC = pd.read_csv(pathW+'Rate2.csv', header=None)
else:
    data1 = pd.read_csv(pathW+'qcd_Ana_rate1.csv', header=None)
    data2 = pd.read_csv(pathW+'meson_Ana_rate1.csv', header=None)
    data3 = pd.read_csv(pathW+'pbrem_Ana_rate1.csv', header=None)
    #dataC = pd.read_csv(pathW+'comb_Ana_rate1.csv', header=None)
    dataC = pd.read_csv(pathW+'comb.csv', header=None)
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

pl.plot(data1[0], data1[1], 'b', label='SHiP (QCD)')
pl.plot(data3[0], data3[1], 'g', label='SHiP (Pbrem)')
pl.plot(data2[0], data2[1], 'r', label='SHiP (Mesons)')
pl.plot(dataC[0], dataC[1], 'black', label='SHiP (Combined)')

#pl.fill_between(data1[0], data1[1], y2=0, color='red', hatch='.', label='SHiP (QCD)')
#pl.fill_between(data3[0], data3[1], y2=0, color='blue', label="SHiP (Brem)") 
#pl.fill_between(data2[0], data2[1], y2=0, color='green', label="SHiP (Mesons)")

#pl.fill_between(data3[0], data3[1], where=data3[0] > y3  , color='green', label="SHiP (Brem)")

pl.fill_between(dataC1[0], dataC1[1], y2=y2, color="grey", label="Excluded Region")
pl.fill_between(dataC2[0], dataC2[1], y2=0, color="grey")
 
pl.xlabel("M [MeV]",fontsize=18)
pl.ylim(1e-10, 1e-2)
pl.xlim(13, 1e4)
pl.ylabel(r"$\epsilon$", fontsize=18)
pl.legend(prop={'size': 10})
pl.show()
