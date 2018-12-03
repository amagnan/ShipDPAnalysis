from array import array
import os,sys,getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], "a:p:d:")

except getopt.GetoptError:
    print 'no file'
    sys.exit()

for o,a in opts:
    if o in ('-a',): analysis = a
    if o in ('-p',): prod = a
    if o in ('-d',): date= a

pathR = "../data/"+date+"/"

f1=open(pathR+prod+'_'+analysis+'_true_table.dat','r')
f2=open(pathR+prod+'_'+analysis+'_e_table.dat','r')
f3=open(pathR+prod+'_'+analysis+'_mu_table.dat','r')
if analysis=='Dau': f4=open(pathR+prod+'_'+analysis+'_tau_table.dat','r')
f5=open(pathR+prod+'_'+analysis+'_pi_table.dat','r')
f6=open(pathR+prod+'_'+analysis+'_ka_table.dat','r')
f7=open(pathR+prod+'_'+analysis+'_4pi_table.dat','r')
f8=open(pathR+prod+'_'+analysis+'_3pi_table.dat','r')
f9=open(pathR+prod+'_'+analysis+'_2pi0_table.dat','r')
f10=open(pathR+prod+'_'+analysis+'_oth_table.dat','r')

f12=open(pathR+prod+'_'+analysis+'_single_table.dat','r')
h=open(pathR+prod+'_'+analysis+'.txt','w')
a1=f1.readlines() 
a2=f2.readlines()
a3=f3.readlines()
if analysis=='Dau': a4=f4.readlines()
a5=f5.readlines()
a6=f6.readlines()
a7=f7.readlines()
a8=f8.readlines()
a9=f9.readlines()
a10=f10.readlines()

a12=f12.readlines()

ev,sig,ve,ang=0.,0.,0.,0.
flag=0
for i,mass in enumerate(a1):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('tot %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])

ev,sig,ve,ang=0.,0.,0.,0.   
flag=0
for i,mass in enumerate(a2):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('e %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])


ev,sig,ve,ang=0.,0.,0.,0. 
flag=0
for i,mass in enumerate(a3):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('mu %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])

if analysis=='Dau':
    ev,sig,ve,ang=0.,0.,0.,0. 
    flag=0
    for i,mass in enumerate(a4):
        k=mass.split(' ')
        if i==0:
            ev+=float(k[2])
            sig+=float(k[3])
            ve+=float(k[4])
            ang+=float(k[5])
            flag=float(k[0])
        elif float(k[0])==flag:
            ev+=float(k[2])
            sig+=float(k[3])
            ve+=float(k[4])
            ang+=float(k[5])
            flag=float(k[0])
        else:
            h.write('tau %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
            h.write('\n') 
            ev,sig,ve,ang=0.,0.,0.,0.
            flag=float(k[0])
            ev+=float(k[2])
            sig+=float(k[3])
            ve+=float(k[4])
            ang+=float(k[5])


ev,sig,ve,ang=0.,0.,0.,0. 
flag=0
for i,mass in enumerate(a5):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('pi %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])


ev,sig,ve,ang=0.,0.,0.,0. 
flag=0
for i,mass in enumerate(a6):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('ka %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])


ev,sig,ve,ang=0.,0.,0.,0. 
flag=0
for i,mass in enumerate(a7):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('4pi %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])


ev,sig,ve,ang=0.,0.,0.,0. 
flag=0
for i,mass in enumerate(a8):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('3pi %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])

ev,sig,ve,ang=0.,0.,0.,0. 
flag=0
for i,mass in enumerate(a9):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('2pi0 %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])

ev,sig,ve,ang=0.,0.,0.,0.
flag=0
for i,mass in enumerate(a10):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('oth %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
 

ev,sig,ve,ang=0.,0.,0.,0.
flag=0
for i,mass in enumerate(a12):
    k=mass.split(' ')
    if i==0:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    elif float(k[0])==flag:
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])
        flag=float(k[0])
    else:
        h.write('single %.8g %.8g %.8g %.8g %.8g'%(flag,ev,sig,ve,ang))
        h.write('\n') 
        ev,sig,ve,ang=0.,0.,0.,0.
        flag=float(k[0])
        ev+=float(k[2])
        sig+=float(k[3])
        ve+=float(k[4])
        ang+=float(k[5])


f1.close() 
f2.close()
f3.close()
if analysis=='Dau': f4.close()
f5.close()
f6.close()
f7.close()
f8.close()
f9.close()
f10.close()

f12.close()
h.close()
