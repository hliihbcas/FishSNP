#!/usr/bin/env python

from lzhanglib.base import Output
import fileinput
import numpy as np
from scipy import stats

def chiqTest(fn,out,co):
	inFile = fileinput.input(fn)
	tt=fn[0].split(".")
	outFileA = Output(".".join([out+"_A","chiqtest"]))
	outFileB = Output(".".join([out+"_B","chiqtest"]))
	outFileD = Output(".".join([out+"_D","chiqtest"]))
	outrawA1 = Output(".".join([out+"_A.1",tt[-1]]))
	outrawA2 = Output(".".join([out+"_A.2",tt[-1]]))
	outrawB3 = Output(".".join([out+"_B3.7",tt[-1]]))
	outrawD1 = Output(".".join([out+"_D1.10",tt[-1]]))
	outrawD2 = Output(".".join([out+"_D2.15",tt[-1]]))
	outSum = Output(".".join([out+"_summary","txt"]))
	co = float(co)
	##Add Title
	outFileB.write("\t".join(["ID","Type","Pvalue(a:2ab:b)","Pvalue(a:b)","Pvalue(a:2ab)","Pvalue(2ab:b)","Counts of each genotypes"])+"\n")
	outFileD.write("\t".join(["ID","Type","Pvalue(1:1)","Counts of each genotypes"])+"\n")
	outFileA.write("\t".join(["ID","Type","Pvalue(1:1:1:1)","Counts of each genotypes"])+"\n")
	outSum.write("\t".join(["Type","Total_subBin","Pass_subBin"])+"\n")
	summary={
		"A.1":["A.1",0,0],
		"A.2":["A.2",0,0],
		"B3.7":["B3.7",0,0],
		"D1.10":["D1.10",0,0],
		"D2.15":["D2.15",0,0]
	}

	##Check input format
	if tt[-1]=="raw":
		flag=2
	elif tt[-1]=="wraw":
		flag=3
	##
	for i in inFile:
		i = i.strip()
		line = i.split("\t")
		genotype = line[flag:]
		gr=line[flag-1]
		genocount = countGenotype(gr,genotype)
		pvalue = chiqtest(gr,genocount)
		summary[gr][1]+=1
		if gr=="B3.7":
			outFileB.write("\t".join([line[0],line[1]]+list(map(str,pvalue))+[":".join(map(str,genocount))+"|"+":".join(["a","ab","b"])])+"\n")
			if (min(pvalue) >= co):
				summary[gr][2]+=1
				outrawB3.write(i)
				outrawB3.write("\n")
		elif gr=="D2.15" or gr=="D1.10":
			outFileD.write("\t".join([line[0],line[1],str(pvalue),":".join(map(str,genocount))+"|"+":".join(["a","ab"])])+"\n")
			if  pvalue >= co:
				summary[gr][2]+=1
				if gr=="D2.15":
					outrawD2.write(i+"\n")
				else:
					outrawD1.write(i+"\n")
		elif gr=="A.1":
			outFileA.write("\t".join([line[0],line[1],str(pvalue),":".join(map(str,genocount))+"|"+":".join(["ac","ad","bc","bd"])])+"\n")
			if pvalue >=co:
				summary[gr][2]+=1
				outrawA1.write(i+"\n")
		elif gr=="A.2":
			outFileA.write("\t".join([line[0],line[1],str(pvalue),":".join(map(str,genocount))+"|"+":".join(["a","ac","ba","bc"])])+"\n")
			if pvalue >= co:
				summary[gr][2]+=1
				outrawA2.write(i+"\n")

	for key in summary.keys():
		outSum.write("\t".join(map(str,summary[key]))+"\n")
	outFileD.close()
	outrawA1.close()
	outrawA2.close()
	outrawB3.close()
	outrawD1.close()
	outrawD2.close()
	outSum.close()
	

def countGenotype(gr,genotype):
	if gr == "B3.7":
		a=genotype.count("a")
		ab=genotype.count("ab")
		b=genotype.count("b")
		return [a,ab,b]
	elif gr == "D1.10" or gr == "D2.15":
		a=genotype.count("a")
		ab=genotype.count("ab")
		return [a,ab]
	elif gr == "A.1":
		ac=genotype.count("ac")
		ad=genotype.count("ad")
		bc=genotype.count("bc")
		bd=genotype.count("bd")
		return [ac,ad,bc,bd]
	elif gr == "A.2":
		a=genotype.count("a")
		ac=genotype.count("ac")
		ba=genotype.count("ba")
		bc=genotype.count("bc")
		return [a,ac,ba,bc]

def chiqtest(gr,one):
	al=sum(one)
	if gr== "B3.7" :
		a=one[0]
		ab=one[1]
		b=one[2]
		ori=al/4
		##a:2b:b
		pvalue1=stats.chi2_contingency(np.array([one,[ori,2*ori,ori]]))[1]
		##a:b
		if a==b and b==0:
			pvalue2=0
		else:
			pvalue2=stats.chi2_contingency(np.array([[a,b],[ori,ori]]))[1]
		##a:ab
		if a==ab and a==0:
			pvalue3=0
		else:
			pvalue3=stats.chi2_contingency(np.array([[a,ab],[ori,2*ori]]))[1]
		##ab:b
		if ab==b and b==0:
			pvalue4=0
		else:
			pvalue4=stats.chi2_contingency(np.array([[ab,b],[2*ori,ori]]))[1]
		##All
		pvalue=[pvalue1,pvalue2,pvalue3,pvalue4]
	elif gr == "D1.10" or gr == "D2.15":
		ori=al/2
		two=[ori,ori]
		f_obs=np.array([one,two])
		pvalue=stats.chi2_contingency(f_obs)[1]
	elif gr == "A.1" or gr == "A.2":
		ori=al/4
		pvalue=stats.chi2_contingency(np.array([one,[ori,ori,ori,ori]]))[1]
	return pvalue


def main(fn,out,co):
	chiqTest(fn,out,co)


if __name__ == '__main__':
		from optparse import OptionParser
		import sys
		usage="\n  %prog [options] [arg1 arg2 ...] inputfile"
		parser=OptionParser(usage,version="%prog [V1.0.1]")
		#parser.add_option("-v","--version", action="store_true", dest="verbose",help="Print version information")
		parser.add_option("-c","--cutoff", dest="cutoff", default=0.01, help="The cutoff value of Chiqtest. Default is 0.01")		
		parser.add_option("-o","--out", dest="out", default='out',help="The prefix of the output files. The default is \"out\"")
		#verbose="/n=====V1.0.0=====/n"
		(options, files) = parser.parse_args()
		if not files:
			parser.print_help()
			print("\n####The input file is needed!####")
			sys.exit(0)
		#start
		main(files,options.out,options.cutoff)
