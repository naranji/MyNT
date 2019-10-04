
'''
Created on Jun 11, 2014

@author: chenani
'''

import os,sys,fnmatch
import pickle as pkl
import numpy as np
import matplotlib.pyplot as pl
import sws_rem_compare as swrem
##########################################
animalPath = "/home/chenani/DATA-clone/G9588"
##########################################FUNCTIONS
def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
        supplied root directory.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield [path,filename]
def save(report, filename):
        '''
        This function saves the ephys object to a .ephys file!
        '''
        pkl.dump(report, open(filename + ".rpt", 'wb'), pkl.HIGHEST_PROTOCOL)
        print "data has been saved to %s.rpt" % filename
##########################################
ephysList = []
vr_ripp_sig = np.array([])
lt_ripp_sig = np.array([])
for ephys in locate("*.ephys",animalPath):
    ephysList.append(ephys)
#########################################
for item in ephysList:
    if str(item).find('leep2') > 1:
        print os.path.join(item[0],item[1])
        ephysFile = os.path.join(item[0],item[1])
        rptFile = item[0]+'/catalog/report.rpt'
        ephys = pkl.load(open(ephysFile, "rb"))
        report = pkl.load(open(rptFile, "rb"))
        csc_rank = report['quality_rank']+report['strength_rank']
        csc = ephys.LFPs[csc_rank.argmin()]
#         csc.filter(140,200)
#         csc.SWS_signal(True)
#         csc.REM_signal()
#         csc.hilbertTransform()
#         csc.ripple_recorder(removeREMripples=True)
#         ephys.save(ephysFile.split('.')[0])
        if str(item).find('/vr/') > 1:
            vr_ripp_sig = np.append(vr_ripp_sig,np.hstack(csc.SWRmix))
            print 'item added to VR'
        else:
            lt_ripp_sig = np.append(lt_ripp_sig,np.hstack(csc.SWRmix))
            print 'item added to real'
x,xfft,y,yfft = swrem.fft(lt_ripp_sig,vr_ripp_sig, csc.sampleFreq, display=True,step=True)
     

        