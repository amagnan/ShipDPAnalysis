from array import array
f=open('excl.dat','r')
l=f.readlines()
s=[]
for x in l:
    a=x.replace(",","")
    y=a.split(" ")
    res=pow(float(y[1]),2.)
    s.append(res)
    print res
