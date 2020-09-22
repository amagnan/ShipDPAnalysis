import os,sys,getopt
import ROOT as r
import numpy as np
from array import array
prod=0
try:
    opts, args = getopt.getopt(sys.argv[1:], "a:p:d:")

except getopt.GetoptError:
    print 'no file'
    sys.exit()

for o,a in opts:
    if o in ('-a',): analysis = a
    if o in ('-p',): prod = a
    if o in ('-d',): date = a

pathR = "../data/"+date+"/"
pathW = "../Exclusion/"
eps,rate=array('d'),array('d')

if prod:
    f=open(pathR+prod+'_'+analysis+'.dat','r')
    k=open(pathR+prod+'_'+analysis+'.csv','w')
 
a=f.readlines()
c1  = r.TCanvas('c1', '',1920,1080)  
gr={}
h=0
Mass=0
Ml,Mu,El,Eu,Rl,Ru=[],[],[],[],[],[]
for i,dummy in enumerate(a):
    j=dummy.split(' ')
    if i==0:
        eps.append(r.TMath.Log10(float(j[1])))
        rate.append(r.TMath.Log10(float(j[2])))
        Mass=float(j[0])
    elif Mass==float(j[0]):
        eps.append(r.TMath.Log10(float(j[1])))
        rate.append(r.TMath.Log10(float(j[2])))
        Mass=float(j[0])
    else:
        gr[h]=r.TGraph(len(eps),eps,rate)
        gr[h].SetTitle("Mass of "+str(Mass)+";Log(#varepsilon);Log(Rate)")
        c1.Update()
        gr[h].SetMarkerStyle(20)
        gr[h].SetMarkerColor(2)
        gr[h].SetLineWidth(3)
        gr[h].Draw("AP")
        ind=rate.index(max(rate))
        #print Mass, eps
        line=r.TLine(max(eps),r.TMath.Log10(2.3),min(eps),r.TMath.Log10(2.3))
        line.SetLineColor(3)
        line.Draw("same")
        #c1.Print(pathR+'pdf/'+prod+'_'+analysis+'_'+str(h)+".pdf")
        ind=rate.index(max(rate))
        for mc in range(-10000000,int(eps[ind]*1000000+1)): 
            mc=mc/1000000.
            rate_exp=gr[h].Eval(mc)
            #if 10**rate_exp<2.34999999999 and 10**rate_exp>2.24999999999:
            if abs(round(10**rate_exp,2)-2.3)<0.07:
                Ml.append(Mass)
                El.append(10**mc)
                Rl.append(abs(10**rate_exp-2.3)) 
        for mc in range(int(eps[ind]*1000000),-2000001):
            mc=mc/1000000.
            rate_exp=gr[h].Eval(mc)
            #if 10**rate_exp<2.34999999999 and 10**rate_exp>2.24999999999:
            if abs(round(10**rate_exp,2)-2.3)<0.07:
                Mu.append(Mass)
                Eu.append(10**mc)
                Ru.append(abs(10**rate_exp-2.3))
 
       
       
        if len(Ru)!=0:
            Iu = Ru.index(min(Ru))
            print 'lowest error', Mass, Eu[Iu]
            k.write('%.13g,  %.13g, %.13g' %(Mass*1000.,Eu[Iu], min(Ru)))
            k.write('\n')
            line2=r.TLine(r.TMath.Log10(Eu[Iu]),-10,r.TMath.Log10(Eu[Iu]),10)
            line2.SetLineColor(3)
            line2.Draw("same")
            bU=np.argsort(Ru) 
            for i,n in enumerate(bU):
                print Mass, Eu[n], Ru[n]
                if i>10: break
 
        if len(Ru)==0 and len(Rl)!=0:
            print 'High eps empty',Mass, max(El)
            print 'High eps empty',Mass, min(El)
        
        if len(Rl)!=0:
            Il = Rl.index(min(Rl))
            print 'lowest error', Mass, El[Il]
            k.write('%.13g,  %.13g, %.13g' %(Mass*1000.,El[Il], min(Rl)))
            k.write('\n')
            line3=r.TLine(r.TMath.Log10(El[Il]),-10,r.TMath.Log10(El[Il]),10)
            line3.SetLineColor(3)
            line3.Draw("same")
            bL=np.argsort(Rl)
            for i,n in enumerate(bL):
                print Mass, El[n], Rl[n]
                if i>10: break
 
        if len(Rl)==0 and len(Ru)!=0: 
            print 'Low eps empty', Mass, max(Eu)
            print 'Low eps empty', Mass, min(Eu)
 
        Ml,Mu,El,Eu,Rl,Ru=[],[],[],[],[],[]
        #c1.Update()
        c1.Print(pathR+'pdf/'+prod+'_'+analysis+'_'+str(h)+".pdf")
        h+=1
        eps,rate=array('d'),array('d')
        eps.append(r.TMath.Log10(float(j[1])))
        rate.append(r.TMath.Log10(float(j[2])))
        Mass=float(j[0])
        #k.close()
