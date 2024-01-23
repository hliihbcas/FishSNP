#!/usr/bin/env python3
import os, sys
# os.environ['MODIN_ENGINE']='dask'
# from distributed import Client
# from collections import Counter
import re
import pandas as pd
# import modin.pandas as pd
import numpy as np
np.set_printoptions(threshold=np.inf)
pd.set_option('display.width',10000)
pd.set_option('max_colwidth',1000)
pd.set_option('display.max_rows', None) # ______________
pd.set_option('display.max_columns', None) # ________________None____________

def vcf (INPUT, OUTPUT, f, s, e, a):
    v_list = []      # ______
    guding_lie = 9   # ____________
    # print(INPUT, OUTPUT, f, s, e, a )
    if INPUT == "-":
        f1 = sys.stdin
    else:
        f1 = open(INPUT)
    f2 = open(OUTPUT, "w")
    lie = 0
    for eachline in f1:
        if eachline[0] == "#" and eachline[1] == "#":  # __________
            lie = lie + 1
            if a:
                f2.writelines(eachline)
        elif eachline[0] =="#" and eachline[1] != "#":  # ________________________index
            #f2.writelines(eachline)
            values = re.split('\s+',eachline)
            # print(values)
            for i in range(0, guding_lie):
                v_list.append(values[i])
            # print(v_list)
            break
    # print(lie)
    data = pd.read_csv(INPUT, sep = '\s+', skiprows=lie-1, chunksize=500000)
    # mylist = []
    # for chunk in data:
    #     mylist.append(chunk)
    # temp_df = pd.concat(mylist, axis=0)
    # del mylist
    # print(data,file=f2)
    if s == "-":
        ll = sys.stdin
    if e == "-":
        l2 = sys.stdin
    if f == "-":
        f3 = sys.stdin
    else:
        f3 = open(f)
    if f:
        samples = []
        for eachline3 in f3:
            split = eachline3.split()
            v_list = v_list + split
        # print(v_list)
        # print(data[v_list], file=f2)
        linshi = 1
        for chunk in data:
            # print(chunk.columns.values)
            # print(chunk.shape)
            # print(chunk)
            if linshi == 1:
                chunk[v_list].to_csv(f2, sep='\t', index=None, line_terminator="\n")
                linshi = linshi + 1
            else:
                chunk[v_list].to_csv(f2, sep='\t', header=0, index=None, line_terminator="\n")
        # temp_df[v_list].to_csv(f2, sep='\t', index = None, line_terminator="\n")
        f2.close()
    elif s:
        if s == "-":
            for eachline3 in ll:
                samples = eachline3.strip().split(",")
                v_list = v_list + samples
                break
        else:
            samples = s.split(",")
            # print(samples)
            v_list = v_list + samples
            # print(v_list)
            # print(data[v_list], file=f2)
            linshi = 1
            for chunk in data:
                if linshi == 1:
                    chunk[v_list].to_csv(f2, sep='\t', index=None, line_terminator="\n")
                    linshi = linshi + 1
                else:
                    chunk[v_list].to_csv(f2, sep='\t', header =0,index=None, line_terminator="\n")
            # temp_df[v_list].to_csv(f2, sep='\t', index=None, line_terminator="\n")
            f2.close()
    elif e:
        if e == "-":
            for eachline3 in l2:
                samples = eachline3.strip().split(",")
                break
            vv = []
            for item in values:
                if item not in samples:
                    vv.append(item)
                    # print(vv)
        else:
            samples = e.split(",")
            # print(samples)
            vv = []
            for item in values:
                if item not in samples:
                    vv.append(item)
            vv.pop()
            # print(vv)
            linshi = 1
            for chunk in data:
                if linshi == 1:
                    chunk[vv].to_csv(f2, sep='\t', index=None, line_terminator="\n")
                    linshi = linshi + 1
                else:
                    chunk[vv].to_csv(f2, sep='\t', header=0, index=None, line_terminator="\n")
            # temp_df[vv].to_csv(f2, sep='\t',index=None, line_terminator="\n")
            f2.close()
    if f != "-":
        f3.close()

def arg(argv):
    import argparse   # ____argparse__
    parser = argparse.ArgumentParser(description="Purpose:\n\
    This script is used to deal with file.vcf. You can sort the sample column by sample name, and select the sample columns which you want.\n\
    ", epilog="For example:\n\
    [1]python3 vcf_sort_filtrate.py test.vcf f.vcf -f sample.txt\n\
    [2]python3 vcf_sort_filtrate.py test.vcf s.vcf -s 1,4,8,98,34,78\n\
    [3]python3 vcf_sort_filtrate.py test.vcf e.vcf -e 1,4,8,98,34,78\n\
    [4]python3 vcf_sort_filtrate.py test.vcf a.vcf -a or python3 vcf_sort_filtrate.py test.vcf e.vcf -e 1,4,8,98,34,78 -a")
    parser.add_argument('INPUT', nargs='?', help="input folder name/path, \"-\" for stdin.")
    parser.add_argument('OUTPUT', type=str, help="output folder name/path")
    parser.add_argument('-f', nargs='?', default=False, type=str, metavar="sample_file",
                        help="select the sample column according to the sample list, \"-\" for stdin. If you want to select the sample columns according to a sample list, please arrange the samples' name into a column or separate the samples' name with ' ' in your sample list. [Default:False]")
    parser.add_argument('-s', nargs='?', default=False, type=str, metavar="sample_line",
                        help="select the sample columns according to the samples' name after the parameter s, \"-\" for stdin. If you want to select the sample columns according to the samples' name after the parameter s, please separate the samples' name with ','. [Default:False]")
    parser.add_argument('-e', nargs='?', default=False, type=str, metavar="sample_line",
                        help="delete the sample columns according to the samples' name after the parameter e, \"-\" for stdin. If you want to delete the sample columns according to the samples' name after the parameter e, please separate the samples' name with ','. [Default:False]")
    parser.add_argument('-a', default=True, action='store_false',
                        help="If you do not want annotation, you can use this parameter. [Default:True]")
    args = parser.parse_args()
    # print(args.INPUT, args.OUTPUT, args.f, args.s, args.e, args.a)
    vcf(args.INPUT, args.OUTPUT, args.f, args.s, args.e, args.a)


if __name__ == '__main__':
    arg(sys.argv)

