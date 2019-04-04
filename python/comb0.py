from array import array
import math
import os,sys,getopt

leptophilic=0 

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:")

except getopt.GetoptError:
    print 'decay to all SM particles'

for o,a in opts:
    if o in ('-d',): date = a

pathW = "../data/"+date+"/"

data1  = open(pathW+'qcd_Ana_rate.dat')
data2  = open(pathW+'meson_Ana_rate.dat')
data3  = open(pathW+'pbrem_Ana_rate.dat')

sum1   = open(pathW+'meson_Ana_sum.dat','r')
sum2   = open(pathW+'pbrem_Ana_sum.dat','r')
sum3   = open(pathW+'qcd_Ana_sum.dat','r')
Rate   = open(pathW+'Rate1.dat','w')

line_1 = data1.readlines()
line_2 = data2.readlines()
line_3 = data3.readlines()
 
line_s1 = sum1.readlines()
line_s2 = sum2.readlines()
line_s3 = sum3.readlines()

def find_Rate(lines,mass,eps):
    for i in lines:
        k = i.replace('\n','')
        #print k
        k = k.split(' ')
        if abs(math.log10(float(k[1])) - math.log10(eps)) <0.1 and abs(float(k[0]) - mass)<0.0001 and float(k[1])!=0.0: return float(k[2])
    else: return 0

def find_N(lines,eps,mass):
    for i in lines:
        k = i.replace('\n','')
        k = k.split(' ')
        if abs(math.log10(float(k[1])) - math.log10(eps)) <0.1 and abs(float(k[0]) - mass)<0.0001: return float(k[4])
    else:
        print "something is wrong"
        return 0

for i in line_2:#pbrem just in case
    i = i.replace('\n','')
    i = i.split(' ')
    if float(i[1])!=0.0: p = find_Rate( line_3, float(i[0]), float(i[1]) )
    if float(i[1])==0.0:p=0
    if p:
        #N3 = find_N( line_s3, float(i[0]),float(i[1]) )
        #N2 = find_N( line_s2, float(i[0]), float(i[1]) )
        #Rate_tot = (float(i[2])*N3 + p*N2)/(N3+N2)
        Rate_tot = float(i[2]) + p
        Rate.write('%.4g %.8g %.8g' %(float(i[0]), float(i[1]), Rate_tot))
        Rate.write('\n')
    #else:
        #Rate_tot = float(i[2])
    #Rate.write('%.4g %.8g %.8g' %(float(i[0]), float(i[1]), Rate_tot))
    #Rate.write('\n')

"""for i in line_3:#qcd rates with pbrem shared
    i = i.replace('\n','')
    i = i.split(' ')
    m = find_Rate( line_2, float(i[0]), float(i[1]) )
    q = find_Rate( line_1, float(i[0]), float(i[1]) )
    if not q and not m:
        Rate_tot = float(i[2])
        Rate.write('%.4g %.8g %.8g' %(float(i[0]), float(i[1]), Rate_tot))
        Rate.write('\n')"""

for i in line_1:#meson rates with pbrem shared
    #Rate_tot=0.
    i = i.replace('\n','')
    #print i
    i = i.split(' ')
    print i[0], i[1]
    if float(i[1])!=0.0: p = find_Rate( line_3, float(i[0]), float(i[1]) )
    if float(i[1])==0.0: p=0
    if p:
        #N1 = find_N( line_s1, float(i[0]),float(i[1]) )
        #N2 = find_N( line_s2, float(i[0]), float(i[1]) )
        #Rate_tot = (float(i[2])*N1 + p*N2)/(N1+N2)
        Rate_tot = float(i[2]) + p
        Rate.write('%.4g %.8g %.8g' %(float(i[0]), float(i[1]), Rate_tot))
        Rate.write('\n')
    #else:
        #Rate_tot = float(i[2])
    #Rate.write('%.4g %.8g %.8g' %(float(i[0]), float(i[1]), Rate_tot))
    #Rate.write('\n')


data1.close()
data2.close()
data3.close()
sum1.close()
sum2.close()
sum3.close()
Rate.close()
