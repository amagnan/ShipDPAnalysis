import pandas as pd
import numpy as np

##### Signal Datasets
data1 = pd.read_csv('Sig.csv', header=None)
#### Constraints
dataC1 = pd.read_csv('Dataset_DP_excluded_upper.csv', header=None)
dataC2 = pd.read_csv('Dataset_DP_excluded_left.csv', header=None)


y2 = np.array([1e-4]*len(dataC1[1]))
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as pl


pl.yscale('log')
pl.xscale('log')
"""for i in range(0,len(data1[0])-1):
    print data1[0][i], data1[1][i]"""

#pl.grid(True)
#pl.grid(alpha=0.1)
pl.fill_between(data1[0], data1[1], y2=0, color='red', label="SHiP")
pl.fill_between(dataC1[0], dataC1[1], y2=y2, color="grey", label="Excluded Region")
pl.fill_between(dataC2[0], dataC2[1], y2=0, color="grey")
pl.xlabel("M [MeV]",fontsize=18)
pl.ylim(1e-20, 1e-4)
pl.xlim(13, 1e4)
pl.ylabel(r"$\epsilon^2$", fontsize=18)
pl.legend(prop={'size': 10})
pl.savefig("Excl.pdf")
pl.show()
