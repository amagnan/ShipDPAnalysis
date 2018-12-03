from array import array
n=raw_input("please enter the name of squared file")
f=open(n,'r')
l=f.readlines()
s=[]
for x in l:
    a=x.replace(",","")
    y=a.split(" ")
    res=pow(float(y[1]),2.)
    s.append(res)
    m=float(y[0])*1000.
    print m, res
