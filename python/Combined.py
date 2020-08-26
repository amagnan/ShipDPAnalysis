from array import array
import math
from decimal import Decimal
import os,sys,getopt

pbremOld=0 
wo=1

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:p:w:")

except getopt.GetoptError:
    print 'decay to all SM particles'

for o,a in opts:
    if o in ('-d',): date = a
    if o in ('-p'): pbremOld = a
    if o in ('-w'): wo = a

pathW = "../data/"+date+"/"

data1M  = open(pathW+'meson_ErrorRateM.dat','w')
data1P  = open(pathW+'meson_ErrorRateP.dat','w')
data3M  = open(pathW+'qcd_ErrorRateM.dat','w')
data3P  = open(pathW+'qcd_ErrorRateP.dat','w')

if not wo:
    data1  = open(pathW+'meson_rate1.dat')
    data3  = open(pathW+'qcd_rate1.dat')
    if not pbremOld:
        data2   = open(pathW+'pbrem_rate1.dat')
        data2M  = open(pathW+'pbrem_ErrorRateM.dat','w')
        data2P  = open(pathW+'pbrem_ErrorRateP.dat','w')
        Rate    = open(pathW+'comb_rate1.dat','w')
        RateM   = open(pathW+'comb_ErrorRateM.dat','w')
        RateP   = open(pathW+'comb_ErrorRateP.dat','w')
    if pbremOld:
        data2   = open(pathW+'pbrem1_rate1.dat')
        data2M  = open(pathW+'pbrem1_ErrorRateM.dat','w')
        data2P  = open(pathW+'pbrem1_ErrorRateP.dat','w')
        Rate    = open(pathW+'comb1_rate1.dat','w')
        RateM   = open(pathW+'comb1_ErrorRateM.dat','w')
        RateP   = open(pathW+'comb1_ErrorRateP.dat','w')

if  wo:
    data1  = open(pathW+'meson_Rate1.dat')
    data3  = open(pathW+'qcd_Rate1.dat')
    if not pbremOld:
        data2   = open(pathW+'pbrem_Rate1.dat')
        data2M  = open(pathW+'pbrem_ErrorRateM.dat','w')
        data2P  = open(pathW+'pbrem_ErrorRateP.dat','w')
        Rate    = open(pathW+'comb_Rate1.dat','w')
        RateM   = open(pathW+'comb_ErrorRateM.dat','w')
        RateP   = open(pathW+'comb_ErrorRateP.dat','w')
    if pbremOld:
        data2   = open(pathW+'pbrem1_Rate1.dat')
        data2M  = open(pathW+'pbrem1_ErrorRateM.dat','w')
        data2P  = open(pathW+'pbrem1_ErrorRateP.dat','w')
        Rate    = open(pathW+'comb1_Rate1.dat','w')
        RateM   = open(pathW+'comb1_ErrorRateM.dat','w')
        RateP   = open(pathW+'comb1_ErrorRateP.dat','w')

line_1 = data1.readlines()#meson
line_2 = data2.readlines()#pbrem
line_3 = data3.readlines()#qcd

err12=0.30
err3P=0.10
err3M=0.20

errC12=(2*err12)**0.5
errM3M=(err12*err12+err3M*err3M)**0.5
errM3P=(err12*err12+err3P*err3P)**0.5

def find_Rate(lines,mass,eps):
    for i in lines:
        k = i.replace('\n','')
        k = k.split(' ')
        if abs(math.log10(Decimal(k[1])) - math.log10(eps)) <0.001 and abs(Decimal(k[0]) - mass)<0.00001:
            return Decimal(k[2])
    else: return 0

for i in line_1:#pbrem just in case
    i = i.replace('\n','')
    i = i.split(' ')
    data1P.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*err12)))
    data1P.write("\n")
    data1M.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*err12)))
    data1M.write("\n")
    p = find_Rate( line_2, Decimal(i[0]), Decimal(i[1]) )
    if p:
        Rate_tot = Decimal(i[2]) + p
        RateP_tot = Rate_tot + Rate_tot*Decimal(errC12)
        RateM_tot = Rate_tot - Rate_tot*Decimal(errC12)
        Rate.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), Rate_tot))
        Rate.write('\n')
        RateP.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), RateP_tot))
        RateP.write('\n')
        RateM.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), RateM_tot))
        RateM.write('\n')

for i in line_2:#pbrem just in case
    i = i.replace('\n','')
    i = i.split(' ')
    data2P.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*err12)))
    data2P.write("\n")
    data2M.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*err12)))
    data2M.write("\n")
    if float(i[0])>0.95 and float(i[0])<1.5:
        Rate.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), Rate_tot))
        Rate.write('\n')
        RateP.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*err12)))
        RateP.write("\n")
        RateM.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*err12)))
        RateM.write("\n")

for i in line_3:#meson rates with pbrem shared
    i = i.replace('\n','')
    i = i.split(' ')
    data3P.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*err3P)))
    data3P.write("\n")
    data3M.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*err3M)))
    data3M.write("\n")
    p = find_Rate( line_2, Decimal(i[0]), Decimal(i[1]) )
    if p:
        Rate_tot = Decimal(i[2]) + p
        RateP_tot = Rate_tot + Rate_tot*Decimal(errM3P)
        RateM_tot = Rate_tot - Rate_tot*Decimal(errM3M)
        Rate.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), Rate_tot))
        Rate.write('\n')
        RateP.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), RateP_tot))
        RateP.write('\n')
        RateM.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), RateM_tot))
        RateM.write('\n')
    if float(i[0])>3.7:
        RateP.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*err3P)))
        RateP.write("\n")
        RateM.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*err3M)))
        RateM.write("\n")


data1.close()
data2.close()
data3.close()
Rate.close()

data1P.close()
data2P.close()
data3P.close()
RateP.close()

data1M.close()
data2M.close()
data3M.close()
RateM.close()
