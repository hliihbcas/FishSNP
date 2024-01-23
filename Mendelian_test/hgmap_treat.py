#!/usr/bin/env python

from lzhanglib.base import Output
import fileinput
import re

def oneMapTreat(fn,out,score,errper,nc):
	inFile = fileinput.input(files=fn)
	outFile = Output(".".join([out,"raw"]))
	#outNAcount = Output(".".join([out,"nacount"]))
	outError= Output(".".join([out+"_error","raw"]))
	outValidsnp= Output(".".join([out+"_validsnp","txt"]))
	outSame = Output(".".join([out+"_same","geno"]))
	# print(outSame)
	outLowc = Output(".".join([out+"_Lowcov","raw"]))
	outParentNA = Output(".".join([out+"_ParentNA","geno"]))
	outInfo = Output(".".join([out+"_info","txt"]))
	outInfo.write("\t".join(["SNP","Type","All_Na","Ori_Na","Ori_Err","score","Ori_Err%"])+"\n")
	outSum = Output(".".join([out+"_summary","txt"]))
	outSum.write("\t".join(["#Type","TotalSNP","ErrSNP","LowcovSNP","ValidSNP","ModifyValidSNP"])+"\n")
	summary = {
		"SAME":["SAME",0],
		"ParentNA":["ParentNA",0],
		"A.1":["A.1",0,0,0,0,0],
		"A.2":["A.2",0,0,0,0,0],
		"D1.9":["D1.9",0,0,0,0,0],
		"D1.10":["D1.10",0,0,0,0,0],
		"D2.14":["D2.14",0,0,0,0,0],
		"D2.15":["D2.15",0,0,0,0,0],
		"B3.7":["B3.7",0,0,0,0,0]
	}
	allSNP=0
	validSNP=0
	lowcovSNP=0
	errSNP=0

	head=inFile.readline()

	for i in inFile:
		i = i.strip()
		marker=markerProp(i)
		
		summary[marker.flag][1]+=1
		# print(summary[marker.flag])
		if marker.flag=="SAME":
			outSame.write(i+"\n")
			outInfo.write("\t".join(map(str,[marker.name,marker.flag,"-","-","-","-","-"]))+"\n")
		elif marker.flag=="ParentNA":
			outParentNA.write(i+"\n")
			outInfo.write("\t".join(map(str,[marker.name,marker.flag,"-","-","-","-","-"]))+"\n")
		else:	
			outInfo.write("\t".join(map(str,[marker.name,marker.flag,int(marker.allnanum),int(marker.nanum),int(marker.errnum),round(marker.score,2),round(marker.errper,2)]))+"\n")
			if marker.errper >= errper :
				outError.write("\t".join(map(str,[marker.name,marker.flag,marker.trans]))+"\n")
				summary[marker.flag][2]+=1
			else:
				if marker.score <= score:
					outLowc.write("\t".join(map(str,[marker.name,marker.flag,marker.trans]))+"\n")
					summary[marker.flag][3]+=1
				else:
					summary[marker.flag][4]+=1
					outValidsnp.write("\t".join(map(str,[marker.name,marker.flag,int(marker.allnanum),int(marker.nanum),int(marker.errnum),round(marker.score,2),round(marker.errper,2)]))+"\n")
					if nc:
						outFile.write("\t".join(map(str,[marker.name,marker.flag2,marker.modify]))+"\n")
					else:
						outFile.write("\t".join(map(str,[marker.name,marker.flag,marker.trans]))+"\n")
	## summary output
	for i in ["A.1","A.2","D1.9","D1.10","D2.14","D2.15","B3.7"]:
		summary[i][5]=summary[i][4]
	if nc:
		summary["D1.9"][5] =0
		summary["D1.10"][5]=summary["D1.10"][4]+summary["D1.9"][4]
		summary["D2.14"][5]=0
		summary["D2.15"][5]=summary["D2.14"][4]+summary["D2.15"][4]

	for i in ["A.1","A.2","D1.9","D1.10","D2.14","D2.15","B3.7"]:
		outSum.write("\t".join(map(str,summary[i]))+"\n")
		validSNP+=summary[i][4]
		lowcovSNP+=summary[i][3]
		errSNP+=summary[i][2]
		allSNP+=summary[i][1]
	allSNP=allSNP+summary["SAME"][1]+summary["ParentNA"][1] 
	outSum.write("\n##########summary##########\n")
	outSum.write(" ".join(["# All SNP number:",str(allSNP)])+"\n")
	outSum.write(" ".join(map(str,["# Homozygous SNP number:",summary["SAME"][1],";",float_percent(summary["SAME"][1],allSNP)]))+"\n")
	outSum.write(" ".join(map(str,["# ParentNA SNP number:",summary["ParentNA"][1],";",float_percent(summary["ParentNA"][1],allSNP)]))+"\n")
	outSum.write(" ".join(map(str,["# Lowcov SNP number:",lowcovSNP,";",float_percent(lowcovSNP,allSNP)]))+"\n")
	outSum.write(" ".join(map(str,["# Err SNP number:",errSNP,";",float_percent(errSNP,allSNP)]))+"\n")
	outSum.write(" ".join(map(str,["# Valid SNP number:",validSNP,";",float_percent(validSNP,allSNP)]))+"\n")
	outSame.close()
	outFile.close()
	outError.close()
	outValidsnp.close()
	outLowc.close()
	outParentNA.close()
	outInfo.close()
	outSum.close()



def float_percent(a,b):
	s=round(100*float(a)/float(b),2)
	return "".join([str(s),"%"])


class markerProp(object):
	def __init__(self,line):
		self.markerInfo=line.split("\t")
		self.name=self.markerInfo[2]
		self.vcfg=self.markerInfo[11:]
		self.allnum=len(self.vcfg)-2
		self.transOri()
		if self.flag!="ParentNA" and self.flag!="SAME":
			self.errCount()
	def transOri(self):
		setP1=set(self.markerInfo[9].split("/"))	
		setP2=set(self.markerInfo[10].split("/"))
		# print(self.markerInfo[9].split("/"))
		# print(self.markerInfo[10].split("/"))
		setPa=setP1|setP2
		numP1=int(len(setP1))
		numP2=int(len(setP2))
		i="\t".join(self.vcfg).replace("/","")
		if ("NA" in setPa):
			self.flag="ParentNA"
			self.trans="\t".join(self.vcfg)
			self.modify=self.trans
		elif(numP1==1 and numP2==1): ##P1 x P2 = aa x bb || aa x aa
			self.flag="SAME"
			self.trans="\t".join(self.vcfg)
			self.modify=self.trans
		else:
			numPa=int(len(setPa))
			maxPa=int(max(list(setPa-set("NA"))))
			if (numPa==2):
				if(numP1==1): ##P1 x P2 = aa x ab D2.15;
					self.flag="D2.15"
					self.flag2="D2.15"
					ir=i.replace(str(list(setP1)[0]),"a").replace(str(list(setP2-setP1)[0]),"b").replace("aa","a").replace("ba","ab").replace("bb","err")
					ir=self.errChange(ir)
					irm=ir
				elif (numP2==1): ##P1 x P2 = ab x aa D1.10;
					self.flag="D1.10"
					self.flag2="D1.10"
					ir=i.replace(str(list(setP2)[0]),"a").replace(str(list(setP1-setP2)[0]),"b").replace("aa","a").replace("ba","ab").replace("bb","err")
					ir=self.errChange(ir)
					irm=ir
				else: ##P1 x P2 = ab x ab B3.7;
					self.flag="B3.7"
					self.flag2="B3.7"
					ir=i.replace(list(setPa)[0],"a").replace(list(setPa)[1],"b").replace("aa","a").replace("ba","ab").replace("bb","b")
					ir=self.errChange(ir)
					irm=ir
			elif (numPa==3):
				if (numP1==1): ##P1 x P2 = cc x ab D2.14; Error Genotyping cc and ab/ba
					self.flag="D2.14"
					self.flag2="D2.15"
					ir=i.replace(str(list(setP1)[0]),"c").replace(str(list(setP2)[0]),"a").replace(str(list(setP2)[1]),"b").replace("ca","ac").replace("cb","bc").replace("cc","err").replace("ba","err").replace("ab","err").replace("bb","err").replace("aa","err")
					ir=self.errChange(ir)
					irm=ir.replace("c","a").replace("aa","a").replace("ba","ab")
				elif (numP2==1): ##P1 x P2 = ab x cc D1.9; Error Genotyping cc and ab/ba
					self.flag="D1.9"
					self.flag2="D1.10"
					ir=i.replace(str(list(setP2)[0]),"c").replace(str(list(setP1)[0]),"a").replace(str(list(setP1)[1]),"b").replace("ca","ac").replace("cb","bc").replace("cc","err").replace("ba","err").replace("ab","err").replace("bb","err").replace("aa","err")
					ir=self.errChange(ir)
					irm=ir.replace("c","a").replace("aa","a").replace("ba","ab")
				else: ##P1 x P2 = ab x ac A.2;
					self.flag="A.2"
					self.flag2="A.2"
					mc=str(list(setP1 & setP2)[0])
					m1=str(list(setP1-(setP1 & setP2))[0])
					m2=str(list(setP2-(setP1 & setP2))[0])
					ir=i.replace(mc,"a").replace(m1,"b").replace(m2,"c").replace("aa","a").replace("ca","ac").replace("ab","ba").replace("cb","bc").replace("bb","err").replace("cc","err")
					ir=self.errChange(ir)
					irm=ir
			elif (numPa==4): ##P1 x P2 = ab x cd A.1; Error Genotyping ab/ba and cd/dc
				self.flag="A.1"
				self.flag2="A.1"
				ir=i.replace(str(list(setP1)[0]),"a").replace(str(list(setP1)[1]),"b").replace(str(list(setP2)[0]),"c").replace(str(list(setP2)[1]),"d").replace("ca","ac").replace("da","ad").replace("cb","bc").replace("db","bd").replace("ba","ab").replace("dc","cd").replace("ab","err").replace("cd","err").replace("bb","err").replace("aa","err").replace("dd","err").replace("cc","err")
				ir=self.errChange(ir)
				irm=ir
			self.trans=ir
			self.modify=irm
	
	def errChange(self,ir):
		irTemp=re.sub("[0-9]\S","err",ir)
		irTemp2=re.sub("\S[0-9]","err",irTemp)
		ir=re.sub("[0-9][0-9]","err",irTemp2)
		return ir

	def errCount(self):
		line=self.trans.split("\t")
		self.nanum=float(line.count("NA"))
		self.errnum=float(line.count("err"))
		self.allnanum=self.nanum+self.errnum
		self.score=1-self.allnanum/self.allnum
		if self.allnum==self.nanum:
			self.errper=0
		else:
			self.errper=self.errnum/(self.allnum-self.nanum)
	

def main(fn,out,score,errper,nc):
	oneMapTreat(fn,out,score,errper,nc)



if __name__ == '__main__':
	from optparse import OptionParser
	import sys
	usage="\n  %prog [options] [arg1 arg2 ...] inputfile"
	parser=OptionParser(usage,version="%prog [V1.0.1]")
	#parser.add_option("-v","--version", action="store_true", dest="verbose",help="Print version information")
	parser.add_option("-s","--score",dest="sc",type="float",default=0.8,help="The genotyping rate of offsprings,0-1,float,default=0.5")
	parser.add_option("-e","--errper",dest="ep",type="float",default=0.05,help="The proportion of false genotyped progenies in successful genotyped offsprings,0-1,float,default=0.05")
	parser.add_option("-n","--nochange", action="store_false", dest="nc",default=True,help="Not change D2.14 into D2.15 and D1.9 into D1.10")
	parser.add_option("-o","--out", dest="out", default='out',help="The prefix of the output files. The default is \"out\"")
	#verbose="/n=====V1.0.0=====/n"
	(options, files) = parser.parse_args()
	if not files:
		parser.print_help()
		print("\n####The input file is needed!####")
		sys.exit(0)
	#start
	main(files,options.out,options.sc,options.ep,options.nc)
