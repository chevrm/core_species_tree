#!/bin/env python

import multiprocessing
import sys, os, re, os.path

def hmmer_mt(fn,bd, pref):
    ## Search
    hmmcmd = 'hmmscan -o '+pref+'.hmmscan.out --tblout '+pref+'.hmmtbl.out -E 1e-5 --cpu 1 --noali '+t+' '+fn
    os.system(hmmcmd)

threads = 10
group = []
for i,a in enumerate(sys.argv):
    if(i==0):
        b = os.path.dirname(os.path.realpath(a))
        continue
    if(i==1):
        t = a
        continue
    pref = os.path.splitext(a)[0]
    pref = pref.split("/")[-1]
    p = multiprocessing.Process(target=hmmer_mt, args=(a,b,pref))
    p.start()
