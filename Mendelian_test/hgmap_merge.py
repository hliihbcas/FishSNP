#!/usr/bin/env python

from lzhanglib.base import Output
import fileinput

def snpMerge(fn,out):
	inFile = fileinput.input(fn)
	outFile = Output(".".join([out+"_snp2subBin","txt"]))
	outRaw = Output("".join([out,"_subBin",".wraw"]))
	merge={}
	## Add Title
	outFile.write("\t".join(["subBinID","subBin_NA_num","subBinWeight","SNPIDs"])+"\n")


	for i in inFile:
		i = i.strip().replace("err","NA")
		line = i.split("\t")
		nanum = str(line.count("NA"))
		content = "\t".join([nanum]+line[1:])
		if content in merge.keys():
			merge[content].append(line[0])
		else:
			merge[content]=[line[0]]
	
	num=1
	for key in merge.keys():
		subBinID="subBinID"+str(num)
		subBinWeight=str(len(merge[key]))
		outFile.write("\t".join([subBinID,key.split("\t")[0],subBinWeight]+merge[key])+"\n")
		outRaw.write("\t".join([subBinID,subBinWeight]+key.split("\t")[1:])+"\n")
		num+=1
	outFile.close()
	outRaw.close()

def main(fn,out):
	snpMerge(fn,out)


if __name__ == '__main__':
	from optparse import OptionParser
	import sys
	usage="\n  %prog [options] [arg1 arg2 ...] inputfile"
	parser=OptionParser(usage,version="%prog [V1.0.1]")
	#parser.add_option("-v","--version", action="store_true", dest="verbose",help="Print version information")
	parser.add_option("-o","--out", dest="out", default='out',help="The prefix of the output files. The default is \"out\"")
	#verbose="/n=====V1.0.0=====/n"
	(options, files) = parser.parse_args()
	if not files:
		parser.print_help()
		print("\n####The input file is needed!####")
		sys.exit(0)
	#start
	main(files,options.out)


