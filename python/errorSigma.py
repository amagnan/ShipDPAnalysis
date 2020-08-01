date='200731'
import math
prods=['comb','comb1','meson','pbrem','pbrem1','qcd']
for prod in prods:
	exec('fP = open("../data/"+date+"/%s"+"_ErrorRateP.dat","w")'%(prod))
	exec('fM = open("../data/"+date+"/%s"+"_ErrorRateM.dat","w")'%(prod))
	exec('f0_%s = open("../data/"+date+"/%s"+"_Rate1.dat","r")'%(prod,prod))
	exec('l0 = f0_%s.readlines() '%(prod))
	for i in l0:
		i=i.replace('\n','')
		i=i.split(' ')
		if "pbrem" in prod:
			fP.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*0.3605551275)))
			fP.write("\n")
			fM.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*0.20)))
			fM.write("\n")
		if "meson" in prod:
			fP.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*0.05)))
			fP.write("\n")
			fM.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*0.05)))
			fM.write("\n")
		if "qcd" in prod:
			fP.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*0.05)))
			fP.write("\n")
			fM.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*0.10)))
			fM.write("\n")
		if "comb" in prod:
			fP.write("%s %s %s"%(i[0],i[1],str(float(i[2])+float(i[2])*0.36742346137194515)))
			fP.write("\n")
			fM.write("%s %s %s"%(i[0],i[1],str(float(i[2])-float(i[2])*0.22912878474779202)))
			fM.write("\n")
	fP.close()
	fM.close()
