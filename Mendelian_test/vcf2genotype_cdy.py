#!/usr/bin/env python3
import os, sys
from collections import Counter
import re
import pandas as pd
# import modin.pandas as pd
import numpy as np

np.set_printoptions(threshold=np.inf)
pd.set_option('display.width', 10000)
pd.set_option('max_colwidth', 1000)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

def vcf(inputfile, outputgeno, s):
    if inputfile == "-":
        f1 = sys.stdin
    else:
        f1 = open(inputfile)
    f2 = open(outputgeno, "w")
    lie = 0
    num = 0
    NAdic = {}
    for eachline1 in f1:
        if eachline1[0] == "#" and eachline1[1] == "#":
            lie = lie + 1
            f2.writelines(eachline1)
        elif eachline1[0] == "#" and eachline1[1] != "#":
            f2.writelines(eachline1)
        else:   #dataframe start
            rows = re.split('\s+',eachline1)
            i = 0
            Chrom1POS = rows[0][:11] + rows[1]
            # print(Chrom1POS)
            na_num = 0
            for row in rows:
                i = i+1
                if i == 1:
                    eachline2 = row + "\t"
                elif i == 2:
                    eachline2 = eachline2 + row + "\t"
                elif i == 3 :
                    eachline2 = eachline2 + Chrom1POS
                elif i>3 and i<10:
                    eachline2 = eachline2 + "\t" + row
                else:
                    if i == 12:
                        keys = row
                    valss = row[:3]
                    if valss == "./.":
                        na_num = na_num + 1
                        num = num + 1
                        eachline2 = eachline2 + "\t" + "NA"
                    elif valss == "0/0" or valss == "0|0":
                        eachline2 = eachline2 + "\t" + "0/0"
                    elif valss == "0/1" or valss == "0|1":
                        eachline2 = eachline2 + "\t" + "0/1"
                    elif valss == "1/0" or valss == "1|0":
                        eachline2 = eachline2 + "\t" + "1/0"
                    elif valss == "1/1" or valss == "1|1":
                        eachline2 = eachline2 + "\t" + "1/1"
            eachline2 = eachline2 + "\n"
            f2.writelines(eachline2)
            NAdic[keys] = na_num
    if s:
        if s == "-":
            f3 = sys.stdin
        else:
            f3 = open(s, "w")
            NAzb = {}
            for k, v in NAdic.items():
                NAzb[k] = v/num*100
                # print(v/num*100)
            print(NAdic, file=f3)
            print(NAzb,file=f3)

def arg(argv):
    import argparse   # 导入argparse包
    parser = argparse.ArgumentParser(description="Purpose:\n\
    This script is used to deal with file.vcf. You can sort the sample column by sample name, and select the sample columns which you want.\n\
    ", epilog="For example:\n\
    [1]python3 vcf2genotype.py test.vcf geno.vcf\n\
    [2]python3 vcf2genotype.py test.vcf l.vcf -s summery_NA.txt. ")
    parser.add_argument('inputfile', nargs='?', help="input folder name/path, \"-\" for stdin.")
    parser.add_argument('outputgeno', type=str, help="output folder name/path")
    parser.add_argument('-s', nargs='?', default=False, type=str, metavar="summary_file",
                        help=" Count the proportion of missing values of samples and markers, and output summary_ NA.txt .Then Draw the frequency distribution and box diagram of Na proportion. [Default:False]")
    args = parser.parse_args()
    # print(args.inputfile, args.outputgeno, args.s)
    vcf(args.inputfile, args.outputgeno, args.s)

if __name__ == '__main__':
    arg(sys.argv)



