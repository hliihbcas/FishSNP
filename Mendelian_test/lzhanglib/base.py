#!/usr/bin/env python
'''
provide basic tools for table-like file operations
'''
import sys
from xqlibs.basTools.versiontk import *  # next, uniopen

def Input(fn, filetype=None, **kw):
    import gzip, bz2, mimetypes
    if fn == '-':
        return sys.stdin
    else:
        if not filetype:
            ## all file types in  /usr/local/lib/python3.8/mimetypes.py line 413-545
            tag=mimetypes.guess_type(fn)
            if (not tag[0] and not tag[1]) or tag[0].startswith('text'): ##(None, None)
                return uniopen(fn, **kw)
            elif not tag[0] and tag[1] == 'gzip': ##(None, 'gzip') gzip file
                return uniopen(gzip.open(fn), **kw)
            elif not tag[0] and tag[1] == 'bzip2': ##(None, 'bzip2') bz2 file
                return uniopen(bz2.open(fn), **kw)
            else:
                print('Err: Unrecognized file type:')
                print(tag)
                sys.exit(1)
                return fn
        else:
            if filetype == 'gzip' or filetype == 'gz':
                return uniopen(gzip.open(fn))
            elif filetype == 'bzip2' or filetype == 'bz2':
                return uniopen(bz2.open(fn))
            elif filetype == 'text' or filetype == 'txt' :
                return uniopen(fn, **kw)
            else:
                print('Err: Unrecognized file type:')
                print(filetype)
                sys.exit(1)
                return fn

def Output(fn, filetype=None, **kw):
    import gzip, bz2
    if not filetype or filetype == 'txt' or filetype == 'text':
        if fn == '-':
            return sys.stdout
        else:
            return uniopen(fn, mode='w', **kw) or sys.stdout
    elif filetype == 'gzip' or filetype == 'gz':
        return uniopen(gzip.open(fn, mode='wb', **kw))
    elif filetype == 'bzip2' or filetype == 'bz2':
        if sys.version < '3':
            return uniopen(bz2.BZ2File(fn, mode='wb', **kw))
        else:
            return uniopen(bz2.open(fn, mode='wb', **kw))
    else:
        print('Error: Unrecognized file type, please ungrade the Output function.')
        sys.exit(1)
        return fn

def sangeropenseq(fn, filetype='abi'):
    from Bio import SeqIO
    a=list(SeqIO.parse(fn, filetype))
    return a[0].seq

def pysamInput(fn, filetype='bam', compressed=False, **kw):
    import pysam
    if isinstance(fn, str) or fn=="-":
        if filetype=='bam' or filetype=='sam':
            return pysam.AlignmentFile(fn, **kw)
            ##**kw  https://pysam.readthedocs.io/en/latest/api.html#sam-bam-cram-files
        elif filetype=='cram':
            return pysam.AlignmentFile(fn, "rc", **kw)
        elif filetype=='fasta' or filetype=='fa' or filetype=='fastq' or filetype=='fq':
            ##**kw  https://pysam.readthedocs.io/en/latest/api.html#fasta-files
            ##**kw  https://pysam.readthedocs.io/en/latest/api.html#fastq-files
            if compressed:
                return pysam.FastxFile(fn, filepath_index_compressed=compressed, **kw)
            else:
                return pysam.FastxFile(fn, **kw)
        elif filetype=='vcf' or filetype=='bcf' or filetype=='vcf.gz':
            ##**kw  https://pysam.readthedocs.io/en/latest/api.html#vcf-bcf-files
            return pysam.VariantFile(fn, **kw)  ##fn can be vcf,bcf,vcf.gz ...
    else:
        print("[ERROR] pysamInput only recognizes filename[STR].")
        sys.exit(1)

def pysamOutput(fn, filetype='bam', **kw):
    import pysam
    ##**kw  https://pysam.readthedocs.io/en/latest/api.html#sam-bam-cram-files
    if isinstance(fn, str) or fn=='-':
        if filetype=='bam':
            return pysam.AlignmentFile(fn, 'wb', **kw)
        elif filetype=='sam':
            return pysam.AlignmentFile(fn, 'w', **kw)
        elif filetype=='cram':
            return pysam.AlignmentFile(fn, 'wc', **kw)
        elif filetype=='vcf' or filetype=='bcf' or filetype=='vcf.gz':
            return pysam.VariantFile(fn, 'w', **kw)
        else:
            return Output(fn, filetype=None, **kw)

def pdLogic(matrix, type=2, logicmethod="and"):  ##input is the pandas dataframe.
    import numpy as np
    import pandas as pd
    if type==1:
        matrix=matrix.T
    if logicmethod=="and":
        func=np.logical_and
    elif logicmethod=="or":
        func=np.logical_or
    elif logicmethod=="not":
        func=np.logical_not
    else:
        print("Wrong logicmethod.")
        sys.exit
    rr=[0]*matrix.shape[0] ##Initiation
    for i in matrix:
        rr=eval("func(rr,np.array(matrix[i]))")
    return rr ##output is a vector of True or False.



#change file into a hash, the first column is the KEY-column.
def Hashin(fn, nl):
    hashi = {}
    with Input(fn) as f_n:
        for i in f_n:
            i = i.strip().split("\t")
            if len(i):
                k = i[nl-1]
                del i[nl-1]
                hashi[k] = i
    return hashi

def Listin(fn):
    with Input(fn) as f_n:
        listA=f_n.read().strip().split('\n')
    return listA


def DictSwap(dictA, duplication=True):
    dictB={}
    if duplication:
        for x, y in dictA.items():
            dictB.setdefault(y, []).append(x)
    else:
        dictB = dict((y, x) for x, y in dictA.items())
    return dictB

def Base2Num(strA, direct=True):
    if direct:
        strA=strA.upper()
        return strA.replace('A', '1').replace('T', '4').replace('C', '2').replace('G', '3')
    else:
        strA=str(strA)
        return strA.replace('1', 'A').replace('4', 'T').replace('2', 'C').replace('3', 'G')

def Pd_bigopen(fn, chunksize=500000, header=None, sep='\t', lineterminator='\n', na_values=['NA', '', ' ']):
    import pandas as pd
    df=[]
    if fn=='-':
        data = sys.stdin if sys.version_info[0] == 2 else sys.stdin.buffer
    else:
        data=fn if isinstance(fn, (int,str,float)) else sys.stdin or fn
    for chunks in pd.read_csv(data, chunksize=chunksize, header=header, sep=sep, lineterminator='\n', na_values=na_values, keep_default_na=False):

        df.append(chunks[1:])
        del chunks
    return df

def PdInput(fn, header=None, seps='\t', linetermin='\n', na_values=['NA', '', ' '], keep_default_na=False, **kw):
    import pandas as pd
    if fn=='-':
        data = sys.stdin if sys.version_info[0] == 2 else sys.stdin.buffer
    else:
        data=fn if isinstance(fn, (int,str,float)) else sys.stdin or fn
    return pd.read_csv(data, header=header, sep=seps, lineterminator=linetermin, na_values=na_values, keep_default_na=keep_default_na, **kw)


def Dd_bighashin(fn, nl):
    import dask.dataframe as dd
    f_n=dd.read_csv(fn)


def Hashlist(fn, st):
    hashi = {}
    with Input(fn) as f_n:
        for i in f_n:
            i = i.strip('\r\n')
            if len(i):
                hashi[i] = st
    return hashi

def Hashmerge(fn, nl):
    hashi = {}
    nl = int(nl)
    with Input(fn) as f_n:
        for i in f_n:
            i = i.strip().split("\t")
            if len(i):
                k = i[nl-1]
                del i[nl-1]
                if(k in hashi):
                    hashi[k]="\t".join([hashi[k],i[0]])  ##one key multiple values
                else:
                    hashi[k]=i[0]
    return hashi

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False




class repeatShow(object):
    '''
    Show a line of string at the same place
    '''
    def __init__(self, msg='', fp=sys.stderr):
        self.fp = fp
        self.backch = '\b'
        if isinstance(msg, (tuple, list)):
            self.msgs = list(msg)
        else:
            self.msgs = [msg]
        self.show()

    def show(self, msg=None, grade=None, append=False):
        'grade can be postive or negative int for location.'
        if msg is None:
            if self.msgs:
                msg = ''.join(self.msgs)
                self.fp.write(msg)
                self.fp.flush()
        mlen = len(self.msgs)
        if grade is not None:
            try:
                grade = int(grade)
            except:
                raise ValueError
        elif append:
            grade = mlen + 1
        else:
            grade = -1
        if grade < 0:
            grade = mlen + grade + 1
        if grade <= 0:
            grade = 1
        if mlen < grade:
            self.msgs.extend(['']*(grade-mlen))
        backs = self.backch * len(''.join(self.msgs[grade-1:]))
        self.msgs[grade-1] = msg
        msg = backs + ''.join(self.msgs[grade-1:])
        self.fp.write(msg)
        self.fp.flush()

    def clear(self):
        if self.msgs:
            self.fp.write(self.backch * len(''.join(self.msgs)))
            self.fp.flush()
            del self.msgs[:]
RepeatShow = repeatShow
