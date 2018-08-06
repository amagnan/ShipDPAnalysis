import pandas as pd
import numpy as np

##### Signal Datasets
data1 = pd.read_csv('Dataset_DP_QCD_theory.csv', header=None)
y2QCD = np.array([2e-16]*len(data1[1]))
data2 = pd.read_csv('Dataset_DP_Mesons2.csv', header=None)
data3 = pd.read_csv('Dataset_DP_Brem2.csv', header=None)
##### Constraints
dataC1 = pd.read_csv('Dataset_DP_excluded_upper.csv', header=None)
dataC2 = pd.read_csv('Dataset_DP_excluded_left.csv', header=None)
dataC3 = pd.read_csv('excl.csv',header=None)

y2 = np.array([1e-4]*len(dataC1[1]))
y3= np.array([1e-20]*len(dataC3[1]))

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as pl
#get_ipython().magic('matplotlib inline')

pl.yscale('log')
pl.xscale('log')

#pl.plot(data1[0], data1[1], label="SHiP Sensitivity")
#pl.fill_between(data1[0], data1[1], y2=y2QCD, where=data1[1]<y2QCD, color='blue', label="SHiP (QCD)")
#pl.fill_between(data1[0], data1[1], y2=y2QCD, where=data1[1]>y2QCD, color='blue', label="SHiP (QCD)")
pl.fill_between(data1[0], data1[1], y2=0,color='blue', label='SHiP (QCD)')
#pl.plot(data1[0], data1[1], linestyle=":", color='blue', label="SHiP (QCD) Theory")
pl.fill_between(data3[0], data3[1], y2=0, color='green', label="SHiP (Brem)")
pl.fill_between(data2[0], data2[1], y2=0, color='red', label="SHiP (Mesons)")
pl.fill_between(dataC1[0], dataC1[1], y2=y2, color="grey", label="Excluded Region")
pl.fill_between(dataC2[0], dataC2[1], y2=0, color="grey")
pl.fill_between(dataC3[0], dataC3[1], y2=y3, color="grey", label="Excluded Region")
#pl.plot(dataC1[0], dataC1[1], label="Excluded Region")
pl.xlabel("M [MeV]",fontsize=18)
pl.ylim(1e-20, 1e-4)
pl.xlim(13, 1e4)
pl.ylabel(r"$\epsilon^2$", fontsize=18)
pl.legend(prop={'size': 10})
pl.show()
