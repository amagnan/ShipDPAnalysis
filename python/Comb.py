from array import array
import math
from decimal import Decimal
import os,sys,getopt

lepto=0 
wo=0

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:l:w:")

except getopt.GetoptError:
    print 'decay to all SM particles'

for o,a in opts:
    if o in ('-d',): date = a
    if o in ('-l'): lepto = a
    if o in ('-w'): wo = a

pathW = "../data/"+date+"/"

if not wo:
    if not lepto:
        data3  = open(pathW+'qcd_rate1.dat')
        data1  = open(pathW+'meson_rate1.dat')
        data2  = open(pathW+'pbrem1_rate1.dat')
        Rate   = open(pathW+'comb1_rate1.dat','w')
    if lepto:
        data3  = open(pathW+'qcd_rate2.dat')
        data1  = open(pathW+'meson_rate2.dat')
        data2  = open(pathW+'pbrem1_rate2.dat')
        Rate   = open(pathW+'comb1_rate2.dat','w')

if wo:
    if not lepto:
        data3  = open(pathW+'qcd_Rate1.dat')
        data1  = open(pathW+'meson_Rate1.dat')
        data2  = open(pathW+'pbrem1_Rate1.dat')
        Rate   = open(pathW+'comb1_Rate1.dat','w')
    if lepto:
        data3  = open(pathW+'qcd_Rate2.dat')
        data1  = open(pathW+'meson_Rate2.dat')
        data2  = open(pathW+'pbrem1_Rate2.dat')
        Rate   = open(pathW+'comb1_Rate2.dat','w')


line_1 = data1.readlines()
line_2 = data2.readlines()
line_3 = data3.readlines()

def find_Rate(lines,mass,eps):
    for i in lines:
        k = i.replace('\n','')
        #print k
        k = k.split(' ')
        if abs(math.log10(Decimal(k[1])) - math.log10(eps)) <0.1 and abs(Decimal(k[0]) - mass)<0.00001: return Decimal(k[2])
    else: return 0

for i in line_1:#pbrem just in case
    i = i.replace('\n','')
    i = i.split(' ')
    p = find_Rate( line_2, Decimal(i[0]), Decimal(i[1]) )
    if p:
        Rate_tot = Decimal(i[2]) + p
        Rate.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), Rate_tot))
        Rate.write('\n')

for i in line_3:#meson rates with pbrem shared
    #Rate_tot=0.
    i = i.replace('\n','')
    #print i
    i = i.split(' ')
    #print i[0], i[1]
    p = find_Rate( line_2, Decimal(i[0]), Decimal(i[1]) )
    if p:
        Rate_tot = Decimal(i[2]) + p
        Rate.write('%.5E %.9E %.9E' %(Decimal(i[0]), Decimal(i[1]), Rate_tot))
        Rate.write('\n')

data1.close()
data2.close()
data3.close()
Rate.close()
