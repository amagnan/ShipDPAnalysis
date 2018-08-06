from array import array
prod=raw_input('qcd,pbrem,meson: ')
f1=open(prod+'_true_table.dat','r')
f2=open(prod+'_e_table.dat','r')
f3=open(prod+'_mu_table.dat','r')
f4=open(prod+'_tau_table.dat','r')
f5=open(prod+'_pi_table.dat','r')
f6=open(prod+'_ka_table.dat','r')
f7=open(prod+'_4pi_table.dat','r')
f8=open(prod+'_3pi_table.dat','r')
f9=open(prod+'_2pi0_table.dat','r')
g=open(prod+'_table.txt','w')
a1=f1.readlines() 
a2=f2.readlines()
a3=f3.readlines()
a4=f4.readlines()
a5=f5.readlines()
a6=f6.readlines()
a7=f7.readlines()
a8=f8.readlines()
a9=f9.readlines()
tot_dp=0
dp,ves,ang=0,0,0
for i in a2:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('e %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
dp,ves,ang=0,0,0
for i in a3:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('mu %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
dp,ves,ang=0,0,0
for i in a4:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('tau %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
dp,ves,ang=0,0,0
for i in a5:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('pi %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
dp,ves,ang=0,0,0
for i in a6:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('ka %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
dp,ves,ang=0,0,0
for i in a7:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('4pi+- %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
dp,ves,ang=0,0,0
for i in a8:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('2pi+-pi0 %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
dp,ves,ang=0,0,0
for i in a9:
    k=i.split(' ')
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('2pi+-2pi0 %i %i %i'%(dp,ves,ang))
g.write('\n')
tot_dp+=dp
n,dp,ves,ang=0,0,0,0
for i in a1:
    k=i.split(' ')
    n+=int(k[2])
    dp+=int(k[3])
    ves+=int(k[4])
    ang+=int(k[5])
g.write('total %i %i %i %i %i'%(n,dp,tot_dp,ves,ang))
g.write('\n')
f1.close() 
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()
f8.close()
f9.close()
g.close()
