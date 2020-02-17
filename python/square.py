from array import array
n=raw_input("please enter the name of squared file: ")
n="../data/200213/"+n+"_FixedRate1.csv"
f=open(n,'r')
l=f.readlines()
s=[]
for x in l:
    x=x.split(",")
    print float(x[1])**2.
