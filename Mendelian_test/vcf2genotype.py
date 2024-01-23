#!/usr/bin/env python
#coding=utf-8

'''
XXXXXX

Staus:
    Being developed ......

Purpose:
    ...

License: GNU General Public License v3.0 (http://www.gnu.org/licenses/gpl-3.0.html)
Author: Dr. Xiao-Qin Xia
Email: xqxia70@gmail.com

Usage:
    XXXXXX  [options]  parameters

Options:
    -h, --help:     Display this message, then exit
    --version:      Print version information
    --profile:      Profile the program to find out the time used by each function
    -i, --input=:   The input file. If not provided, pop the first value in
                    {parameters}, or use STDIN as input. '-' means STDIN
    -o, --output=:  The output file. If not provided, pop the first value in
                    {parameters}, or use STDOUT as output. '-' means STDOUT
    -f, --fish-na=: The output file for NA number of each fish sample
    -p, --pos-na=:  The output file for NA number of each variation position

parameters:
    Input file, Output file in order, in case that some of these augments are not provided as options.

'''

import os, sys, re
from basetools import parseOpts, exitMsg, doc_str, TableFile, InputFile, OutputFile

__doc__ = eval(doc_str)()
__version__ = '0.0.1'

def Main(fsrc, fobj, na_samp, na_pos,  col_samp=10, na='./.', sep='\t', linesep=os.linesep):
    fpsrc = InputFile(fsrc)
    for line in fpsrc:
        if line.startswith('##'):
            continue
        elif not line.startswith('#'):
            exitMsg('No title row in the input file! exit ...')
        else:
            break
    col_samp = col_samp - 1
    titles = line.rstrip('\r\n').split('\t')
    samp_names = titles[col_samp:]
    n_samp = len(samp_names)
    samp_na = [0]*n_samp
    idx_samp = list(range(n_samp))
    fpobj = OutputFile(fobj)
    fpobj.write(line.encode())
    fpos = OutputFile(na_pos)
   # fna=open("test.txt","w")
    fpsrc = TableFile(fpsrc)
    for items in fpsrc:
        samps = list(map(lambda s:s.split(':')[0].replace('|', '/').replace(na,'NA'), items[col_samp:]))
        items[col_samp:] = samps
        n_na = samps.count('NA')
        if n_na:
            list(map(lambda i,s=samps,n=samp_na:s[i]==na and n.__setitem__(i, n[i]+1), idx_samp))
        items[2]=":".join(items[0:2])
        fpos.write(('%s:%s\t%s\n' % (items[0], items[1], n_na)).encode())
        #fpos.write('%s:%s\n' % (items[0].encode(), items[1].encode()))
        #fna.write('%s:%s\n' % (items[0].encode(), items[1].encode()))
        #fna.close
        fpobj.write(('\t'.join(items)+'\n').encode())
    OutputFile(na_samp).write(('\n'.join(map('\t'.join, zip(samp_names, map(str, samp_na))))+'\n').encode())

if __name__ == '__main__':
    from getopt import getopt
    optlst, args = getopt(sys.argv[1:], 'hi:o:f:p:', ['help', 'version', 'input=', 'output=', 'fish-na=', 'pos-na=', 'profile'])

    # print help information
    optdic = dict(optlst)
    if '-h' in optdic or '--help' in optdic:
        exitMsg(__doc__)
    if '--version' in optdic:
        exitMsg(sys.argv[0] + ': ' + __version__ + os.linesep)  # use sys.argv[0] instead of __file__ due to py2exe
    no_profile = '--profile' not in optdic

    err = []
    # options that can have multiple values, e.g., '--rep-opt'
    repopts = [] 
    # mandatory options that musted be supplied, e.g., ['--mand-opt', ('-c', '--col'), ...]
    mandopts = [('-f', '--fish-na'), ('-p', '--pos-na')] 
    # options that will be automatically dealt with. e.g. [(int, ('-n', '--num'), 'num', 3), ...] -- (type, (option names), target variable names, default values). type can be int, float, str, bool, or any self-defined functions, e.g. lambda s:set(re.split(r'[,;|]', s)).
    autopts = [(str, ('-f', '--fish-na'), 'na_samp', ''), (str, ('-p', '--pos-na'), 'na_pos', '')] 
    # tuples for options (in repopts) that should have the same length, e.g., ('gene_chr', 'gene_strand', 'gene_file')
    matchvars = [] 
    # variables with fixed location, will be passed to Main using *fixopts, the rest will be passed to Main using **dictopts
    fixvars = [] # positional variables, maybe this is not necessary.

    optlist, optdict = parseOpts(optlst, repopts=repopts, mandopts=mandopts, autopts=autopts, matchopts=matchvars, fixopts=fixvars, err=err)

    # specific options can be added over here
    try:
        fsrc = optdic.get('-i', optdic.get('--input', None)) or args.pop(0)
        if fsrc == '-':
            fsrc = sys.stdin
    except: 
        # disable this if clause if you are waiting for user's input from TTY
        if os.isatty(file.fileno(sys.stdin)): 
            exitMsg(__doc__)
        fsrc = sys.stdin
    try:
        fobj = optdic.get('-o', optdic.get('--output', None)) or args.pop(0)
        if fobj == '-':
            fobj = sys.stdout
    except: fobj = sys.stdout

    if args:
        err.append('Error - unrecognized parameters:%s\t%s' % (os.linesep, ', '.join(args)))
    # quit on error
    if err:
        err.append('%s%sPlease type "%s -h " for help.%s' % (os.linesep, os.linesep, sys.argv[0], os.linesep))
        exitMsg(err, out=sys.stderr)

    # start job
    if no_profile:
        Main(fsrc, fobj, *optlist, **optdict)
    else:
        import profile, pstats, tempfile
        fh, fn = tempfile.mkstemp(prefix='prof_', suffix='.out')
        #fp = os.fdopen(fh, 'wt')
        profile.run("Main(fsrc, fobj, *optlist, **optdict)", fn)
        p = pstats.Stats(fn)
        p.strip_dirs().sort_stats("cumulative").print_stats()
        os.close(fh)
    
