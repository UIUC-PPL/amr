#!/usr/bin/env python

import sys, os, re, shutil, subprocess

nIterations = int(sys.argv[1])
outFiles = os.listdir('out')

for inum in range(nIterations):
    r = re.compile(r'out_[01]*_'+str(inum)+'$')
    fnames = filter(r.match, outFiles)
    dest = open('out/' + str(inum), 'w')
    for f in fnames:
        fhandle = open('out/' + f, 'r')
        dest.write(fhandle.read())
    dest.flush()
    dest.close()

for inum in range(nIterations):
    subprocess.call(['scp', 'out/'+str(inum), 'charm.cs.uiuc.edu:/www/'])
