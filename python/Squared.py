from array import array
n=raw_input("please enter the name of squared file: ")
f=open(n,'r')
l=f.readlines()
s=[]
for x in l:
    #a=x.replace(""," ")
    x=x.split(" ")
    #res=pow(float(y[1]),1/2.)
    r = float(x[2])*40.6
    #s.append(res)
    #m=float(x[0])
    print x[0], x[1], r
    #print res
