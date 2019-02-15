from array import array
import os,sys,getopt
lept= False
try:
    opts, args = getopt.getopt(sys.argv[1:], "d:p:l", ["date=","production="]) 
except getopt.GetoptError:
    sys.exit()
for o,a in opts:
    if o in ('-l','--leptophilic',): lept = a
    if o in ('-d','--date',): date = a
    if o in ('-p', '--production',): prod = a 
#n=raw_input("please enter the production: ")
#d=raw_input("please enter the date: ")
inp="../data/"+date+"/"+prod+"_Ana_sum.dat"
f=open(inp,'r')
l=f.readlines()
m=[]
e=[]
for x in l:
    y=x.split(" ")
    if int(y[3])<1000 or (float(y[3])-float(y[4]))>= float(y[3])/5.: 
        m.append(y[0])
        e.append(y[1])
        #print y[0], y[1]
    else:
        m.append(0)
        e.append(0)
if not lept: inp2="../data/"+date+"/"+prod+"_Ana_rate1.dat"
else: inp2="../data/"+date+"/"+prod+"_Ana_rate2.dat"
g=open(inp2,'r') 

k=g.readlines()
print 'new rate file'
i=0
for x in k:
    x=x.replace("\n","")
    y=x.split(" ")
    i+=1
    if not (y[2]==0. or (y[0]==m[i-1] and y[1]==e[i-1])): print y[0], y[1], y[2]
