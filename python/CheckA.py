date1='190913'
date2='190516'

prods=['meson','pbrem','qcd']
for prod in prods:
    exec('f1_%s = open("../data/"+date1+"/%s_Rate1.dat","r")'%(prod,prod))
    exec('f2_%s = open("../data/"+date1+"/%s_Rate2.dat","r")'%(prod,prod))
    exec('l1_%s = f1_%s.readlines() '%(prod,prod))
    exec('l2_%s = f2_%s.readlines() '%(prod,prod))
    exec('f01_%s = open("../data/"+date2+"/%s_Rate1.dat","r")'%(prod,prod))
    exec('f02_%s = open("../data/"+date2+"/%s_Rate2.dat","r")'%(prod,prod))
    exec('l01_%s = f01_%s.readlines() '%(prod,prod))
    exec('l02_%s = f02_%s.readlines() '%(prod,prod))

def find(l0,l_i):
    for i in l_i:
        i=i.replace('\n','')
        i = i.split(' ')
        #print i[0],i[1]
        if float(l0[0])==float(i[0]) and float(l0[1])==float(i[1]):
            return True
    return False

def loop(l0,l):
    for i in l0:
        i=i.replace('\n','')
        i=i.split(' ')
        if not find(i,l): print float(i[0]), float(i[1])
        #print i[0], i[1]
 
for prod in prods:
    print prod
    exec('loop(l01_%s,l1_%s)'%(prod,prod))
    #exec('loop(l02_%s,l2_%s)'%(prod,prod))
