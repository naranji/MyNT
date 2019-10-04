'''
Created on Oct 3, 2013

@author: Alireza
'''
import sys
import scipy as sp
# sys.argv[1:] = raw_input('Enter command line arguments: ').split()
# for arg in sys.argv: 
#     print arg
# print sys.argv[0]
# #########Testing of loading a t file
# spikes = []
# header = ''
# with open(sys.argv[1], 'rb') as f:
#      headerFinished =  False
#      while not headerFinished:
#           line = f.readline()
#           if line.startswith('%%ENDHEADER'):
#               headerFinished = True
#           header += line
# header =  header.split('\n')
# for item in header:
#     print item, '\n'
#     
def lpass(x):
    cutoff = 500.
    fs = 44100.
    nyq = fs / 2.
    filterorder = 5
    b, a = sp.signal.filter_design.butter(filterorder, cutoff / nyq)
    return sp.signal.lfilter(b, a, x)
print sys.argv[1]