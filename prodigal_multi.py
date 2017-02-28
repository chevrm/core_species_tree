#!/bin/env python

import multiprocessing
import sys, os, re, os.path

def prodigal_mt(fn,bd, pref):
    ## Train
    traincmd = 'prodigal -t '+pref+'.ptrain -c -i '+fn
    os.system(traincmd)
    ## Gene call
    gccmd = 'prodigal -c -i '+fn+' -a '+pref+'.faa -d '+pref+'.orf -t '+pref+'.ptrain -f gff -o '+pref+'.gff'
    os.system(gccmd)

threads = 10
group = []
for i,a in enumerate(sys.argv):
    if(i==0):
        b = os.path.dirname(os.path.realpath(a))
        continue
    pref = os.path.splitext(a)[0]
    p = multiprocessing.Process(target=prodigal_mt, args=(a,b,pref))
    p.start()
