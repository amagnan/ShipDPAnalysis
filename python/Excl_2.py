import pandas as pd
import numpy as np
import os,sys,getopt
from array import array

leptophilic = 0
comp = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], "l:c:", ["leptophilic=","comp="])

except getopt.GetoptError:
    print 'exit the system'
    sys.exit()

for o,a in opts:
    if o in ('-l',): leptophilic = a
    if o in ('-c',): comp = a

pathR = "../Exclusion/"

if comp==0:
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

elif comp=="pbrem":
    if leptophilic:
        data1 = pd.read_csv(pathR+'qcd_Rate2.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate2.csv', header=None)
        data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate2.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate2.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate2.csv', header=None)
    else:
        data1 = pd.read_csv(pathR+'qcd_Rate1.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate1.csv', header=None)
        data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate1.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate1.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate1.csv', header=None)

elif comp=="meson":
    if leptophilic:
        data1 = pd.read_csv(pathR+'meson1_Rate2.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate2.csv', header=None)
    else:
        data1 = pd.read_csv(pathR+'meson1_Rate1.csv', header=None)
        data2 = pd.read_csv(pathR+'meson_Rate1.csv', header=None)

elif comp=="pbremOnly":
    if leptophilic:
        data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate2.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate2.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate2.csv', header=None)
    else:
        data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate1.csv', header=None)
        dataC = pd.read_csv(pathR+'comb_Rate1.csv', header=None)
        dataO = pd.read_csv(pathR+'comb1_Rate1.csv', header=None)

elif comp=="pbrem0":
    if leptophilic:
        data3 = pd.read_csv(pathR+'pbrem_Rate2.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate2.csv', header=None)
    else:
        data3 = pd.read_csv(pathR+'pbrem_Rate1.csv', header=None)
        data4 = pd.read_csv(pathR+'pbrem1_Rate1.csv', header=None)

else: print comb,'something wrong!!'

##### Constraints
dataC1 = pd.read_csv(pathR+'Dataset_DP_excluded_upper.csv', header=None)
dataC2 = pd.read_csv(pathR+'Dataset_DP_excluded_left.csv', header=None)

y2 = np.array([1e-4]*len(dataC1[1]))

#import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as pl

pl.yscale('log')
pl.xscale('log')

if comp==0:
    pl.plot(data2[0], data2[1], 'red', linewidth=1.8,label='Secondary Meson Decays')
    pl.plot(data3[0], data3[1], 'green', linewidth=1.8,label='Proton Bremsstrahlung')
    pl.plot(data1[0], data1[1], 'blue', linewidth=1.8,label='Quantum ChromoDynamics')
    pl.plot(dataC[0], dataC[1], 'black', linewidth=1.6, label='Combined Production')

if comp=="pbrem":
    pl.plot(data2[0], data2[1], 'red',      linewidth=1.8, label = 'Secondary Meson Decays')
    pl.plot(data4[0], data4[1], 'yellow',   linewidth=1.8, label = 'Proton Bremsstrahlung-noFF')
    pl.plot(data3[0], data3[1], 'green',    linewidth=1.8, label = 'Proton Bremsstrahlung-FF')
    pl.plot(data1[0], data1[1], 'blue',     linewidth=1.8, label = 'Quantum ChromoDynamics')
    pl.plot(dataO[0], dataO[1], 'purple',   linewidth=1.6, label = 'Combined Production-noFF')
    pl.plot(dataC[0], dataC[1], 'black',    linewidth=1.6, label = 'Combined Production-FF')

if comp=='meson':
    pl.plot(data2[0], data2[1], 'red', linewidth=1.8,label='Secondary Meson Decays')
    pl.plot(data1[0], data1[1], 'yellow', linewidth=1.8,label='OLD(Secondary Meson Decays)')

if comp=="pbremOnly":
    pl.plot(data4[0], data4[1], 'yellow',   linewidth=1.8, label = 'Proton Bremsstrahlung-noFF')
    pl.plot(data3[0], data3[1], 'green',    linewidth=1.8, label = 'Proton Bremsstrahlung-FF')
    pl.plot(dataO[0], dataO[1], 'purple',   linewidth=1.6, label = 'Combined Production-noFF')
    pl.plot(dataC[0], dataC[1], 'black',    linewidth=1.6, label = 'Combined Production-FF')

if comp=="pbrem0":
    pl.plot(data4[0], data4[1], 'yellow',   linewidth=1.8, label = 'Proton Bremsstrahlung-noFF')
    pl.plot(data3[0], data3[1], 'green',    linewidth=1.8, label = 'Proton Bremsstrahlung-FF')

pl.fill_between(dataC1[0], dataC1[1], y2=y2, color="grey", label="Excluded Region")
pl.fill_between(dataC2[0], dataC2[1], y2=0, color="grey")

if leptophilic:
    pl.ylim(1.8e-18, 1.1e-7)
    pl.xlim(12, 5.1e3)
if not leptophilic:
    pl.ylim(7e-18, 1e-7)
    pl.xlim(12, 4.1e3)

pl.xlabel("M [MeV]",fontsize=20)
pl.ylabel(r"$\epsilon^2$", fontsize=20)
pl.legend(prop={'size': 16})
if comp: pl.savefig(pathR+"Excl_"+str(comp)+'_'+str(leptophilic)+".pdf")
if not comp:  pl.savefig(pathR+"Excl_"+str(leptophilic)+".pdf")
pl.show()
