'''
Created on 26 Mar 2014

@author: chenani
'''
from __future__ import division
import os, sys
import numpy as np
import matplotlib as mpl
mpl.matplotlib_fname()
import matplotlib.pyplot as pl
import pickle as pkl
# add additional custom paths
extraPaths = ["/home/chenani/pypacks/lib/python2.7/site-packages", \
              "/home/thurley/data/", \
              "/home/haas/packages/lib/python2.6/site-packages", \
    os.path.join(os.path.abspath(os.path.dirname(__file__)), '../scripts')]
for p in extraPaths:
    if not sys.path.count(p):
        sys.path.insert(1, p)
# Additional Modules
import signale, trajectory, Recordings
if len(sys.argv) > 1:
    try: 
        fileName = sys.argv[1]
    except:
        print "Usage:", sys.argv[0], "foo.ephys"; sys.exit(1)
else:
    fileName = raw_input('Enter command line arguments: ').split()[0]
rec = pkl.load(open(fileName, "rb"))

