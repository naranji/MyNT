#! /usr/bin/env python

'''
This Script determines the periods of slow wave sleep(SWS) from LFP and trajectory
author: achenani July 2016
'''

import sys
sys.path.append('/home/chenani/ownCloud/Workspaces/Eclipse/dataAnalysis/Sleep-current/src/')
import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import scipy.signal as scsig
import spectrum as sp
import signale
import cPickle as pkl
import time
import statsmodels.api as sm
import colormaps as cm
import seaborn as sns
import signale.tools as tools
import os
import trajectory
from sklearn.cluster import KMeans
from matplotlib.mlab import specgram
import matplotlib.mlab as mlab


##########FUNCTIONS

def whitenARMA(sig,AR=2,MA=0):
        arma = sm.tsa.ARMA(sig, (AR,MA)).fit(disp=0)
        print 'ARMA parameters calculated for order(%s,%s)' %(AR,MA)
        return arma.resid
def MA(array,ord=2):
    ma = np.array([array[ii-ord/2:ii+ord/2].mean() for ii in range(array.size) if ii >= ord ],dtype=float)
    ma = np.insert(ma,0,array[0:ord/2])
    ma = np.append(ma,array[-(ord/2)-1:-1])
    return ma
def zScore(array):
    return(array - np.average(array))/ np.std(array)
def nvt_loader(filename):
    """ 
    Memory map the Neuralynx .nvt format
    Fields
    -------
    swstx
    swid
    sw_data_size
    qTimeStamps       Cheetah timestamp for this record. This value is in microseconds.
    dwpoints          Points with the color bitfield values for this record.This is a 400
                      element array.  See Video Tracker in reference. 
    sncrc
    dnextracted_x     Extracted X location of the target being tracked. 
    dnextracted_y     Extracted Y location of the target being tracked.
    dnextracted_angle The calculated head angle in degrees clockwise from the positive Y
                      axis. Zero will be assigned if angle tracking is disabled.
    dnTargets         Colored targets using the samebitfield format used by the dwPoints array.
                      Instead of transitions, the bitfield indicates the colors that make up 
                      each particular target and the center point of that target.  This is a 50
                      element array sorted by size from largest (index 0) to smallest(index 49).
                      A target value of 0 means that no target is present in thatindex location.
                      See Video Tracker Bitfield Information in reference. 
    
    Reference:
    ----------
        http://neuralynx.com/software/NeuralynxDataFileFormats.pdf
    """ 
    nev_dtype = np.dtype([
        ('swstx'              , '<i2'),
        ('swid'               , '<i2'),
        ('sw_data_size'       , '<i2'),
        ('qTimeStamps'        , '<u8'),
        ('dwPoints'           , '<u4',(400,)),
        ('sncrc'              , '<i2'),
        ('dnextracted_x'      , '<i4'),
        ('dnextracted_y'      , '<i4'),
        ('dnextracted_angle'  , '<i4'),
        ('dnTargets'          , '<i4',(50,)),
    ])
    return np.memmap(filename, dtype=nev_dtype, mode='readwrite',
       offset=(16 * 2**10))

def mySpecgram(sig,minFreq = None,maxFreq = None,overlap=None,
               windowSize= None,NFFT = None,Fs = None,method=None):
    if not overlap:
        overlap = 48
    kIdx = 0
    timeSlices = []
    timeCenters = []
    while kIdx < sig.size:
        #timeEdges = np.arange(0,Nbins*windowSize,windowSize)
        timeSlices.append([kIdx,kIdx+windowSize])
        kIdx += windowSize - overlap
    Pxx = []
    #print timeSlices
    for item in timeSlices:
        timeCenters.append((item[0]+np.diff(item)/2)[0])
        data = sig[item[0]:item[1]]
        if data.size < windowSize:
            ddd = windowSize - data.size
            if ddd%2 ==0:
                data = np.lib.pad(data, (ddd/2,ddd/2), 'constant', constant_values=(0,0))
            else:
                data = np.lib.pad(data, (ddd/2,ddd/2+1), 'constant', constant_values=(0,0))
        if method == 'yule':
            p = sp.pyule(data, 550, norm='biased', NFFT=NFFT,sampling=Fs)
        elif method =='fft':
            p = sp.Periodogram(data,NFFT=NFFT,sampling=Fs)
        p();
        freqz = np.array(p.frequencies())
        Pxx.append(p.psd)
    Pxx = np.array(Pxx).T
    if not minFreq:
        minFreq = freqz.min()
    if not maxFreq:
        maxFreq = freqz.max()
    indexstart = np.where(freqz >= minFreq)[0][0]
    indexend   = np.where(freqz <= maxFreq)[0][-1]
    freqz = freqz[indexstart:indexend]
    Pxx = Pxx[indexstart:indexend]
    timeCenters = np.array(timeCenters)
    return timeCenters,freqz,Pxx


###########PATHS
#nvtFile = animalPath+'VT1.nvt'
#localPath = '/mnt/Data/ephysdata/Rats/'

try:
    animalPath = sys.argv[1]
    selectionPattern = sys.argv[2]

except IndexError:
    print 'usage: %s path-to-data-directory selection-pattern' %sys.argv[0]
pathSplit = np.array(animalPath.split('/'))
animalidx = np.where(['Rat' in item for item in pathSplit])[0]
animalID  = pathSplit[animalidx][0]
print animalID
###LFP path list
lfpPaths = []
for lfp in tools.locate(selectionPattern,animalPath):
    lfpPaths.append(os.path.join(lfp[0],lfp[1]))
lfpPaths = sorted(lfpPaths)
for path in lfpPaths:
    lfp = pkl.load(open(path,'rb'))
    sessionID   = lfp.tags['file'].split('/')[-2].split('-')[-1]
    sessionName = lfp.tags['file'].split('/')[-3]
    sessionDate = lfp.tags['file'].split('/')[-3].split('_')[-2]
    print '%s is loaded for processing!!!' %path
    nvtFile = "/".join(path.split('/')[:-1]) + '/VT1.nvt'
    ######Load Trajectory
    print 'Loading trajectory file %s' %nvtFile
    ts = np.array(nvt_loader(nvtFile)['qTimeStamps'] / 1e6)
    x = np.array(nvt_loader(nvtFile)['dnextracted_x'])
    y = np.array(nvt_loader(nvtFile)['dnextracted_y'])
    HD_angle = np.array(nvt_loader(nvtFile)['dnextracted_angle'])
    r = np.sqrt(np.power(x, 2) + np.power(y, 2))
    #######calculating Moving average
    rma = MA(r, ord=50)
    ama = MA(HD_angle, ord=50)
    #######Extracting Immobile periods
    cond1 = pd.rolling_std(rma, 50) < 0.7
    cond2 = pd.rolling_std(ama, 50) < 1.5
    for cond in [cond1, cond2]:
        chIdx = np.where(pd.rolling_std(cond, 2) > 0)[0]
        kk = 0
        while kk < chIdx.size:
            initIdx = chIdx[kk - 1]
            finIdx = chIdx[kk]
            if finIdx - initIdx < 300 and finIdx - initIdx > 0:
                cond[initIdx:finIdx] = np.ones(finIdx - initIdx)
            kk += 2
    cnd = np.logical_and(cond1, cond2)
    ####Extracting long spans of immobility!
    sleepList = []
    chIdx = np.where(pd.rolling_std(cnd, 2) > 0)[0]
    if cnd[chIdx[0] - 1]:
        if chIdx[0] > 150:
            sleepList.append(ts[0], ts[chIdx[0]])
            kk = 1
    else:
        kk = 0
    while kk < len(chIdx):
        try:
            chIdx[kk + 1] - chIdx[kk] > 150
            sleepList.append([ts[chIdx[kk]], ts[chIdx[kk + 1]]])
            kk += 2
        except IndexError:
            print kk
            kk += 2

            sleepList = np.array(sleepList)
    ######saving time limits of immobility periods
    ###### Creating Spectrogram of Immobility periods!
    SleepTimeIndex = []
    Pxx = np.zeros([39,1])
    for item in sleepList[:]:
        ii, ti = signale.tools.findNearest(lfp.timeAxis, item[0] * 1e3)
        ff, tf = signale.tools.findNearest(lfp.timeAxis, item[1] * 1e3)
        immoChunkLFP = lfp.signal[ii:ff]
        if immoChunkLFP.size > 2e4:
            t, frq, p = mySpecgram(immoChunkLFP, minFreq=0.1, maxFreq=20, windowSize=2048,
                                   method='yule', Fs=2000, overlap=1024, NFFT=4096)
        SleepTimeIndex.append([ti, tf - ti, t.size])
        Pxx = np.append(Pxx, p, axis=1)
    Pxx = np.delete(Pxx, 0, 1)  # delete the first column of zeros!
    SleepTimeIndex = pd.DataFrame(SleepTimeIndex, columns=['initTime', 'duration', 'durationSeconds'])
    ##### Averaging delta, theta and low gamma bands
    Dmean = frq[frq < 5].mean()
    Tmean = frq[np.logical_and(frq > 6, frq < 12)].mean()
    Gmean = frq[frq > 12].mean()
    Dbounds = np.where(frq < 5)[0][[0, -1]]
    Tbounds = np.where(np.logical_and(frq > 6, frq < 12))[0][[0, -1]]
    Gbounds = np.where(frq > 12)[0][[0, -1]]
    DavgPwr = Pxx[Dbounds[0]:Dbounds[-1] + 1, :].mean(axis=0)
    TavgPwr = Pxx[Tbounds[0]:Tbounds[-1] + 1, :].mean(axis=0)
    GavgPwr = Pxx[Gbounds[0]:Gbounds[-1] + 1, :].mean(axis=0)
    DpwrNormal = DavgPwr * (Dmean / Gmean)
    TpwrNormal = TavgPwr * (Tmean / Gmean)
    ######Clustering sleep states using avg. power in these three different bands.
    pwrData = np.array([DpwrNormal, TpwrNormal, GavgPwr]).T
    Klusters = KMeans(n_clusters=2).fit_predict(pwrData)
    k2 = Klusters.copy()
    ###################################################################
    #####Getting rid if islands up to size 4!!!
    UpIdx = np.where(Klusters == 1)[0]
    DownIdx = np.where(Klusters == 0)[0]

    if UpIdx.size > DownIdx.size:
        SinularIdx = np.where(np.diff(DownIdx) == 2)[0]
        SinularIdx = np.int32(DownIdx[SinularIdx] + np.ones(SinularIdx.size))
        np.put(Klusters, SinularIdx, np.zeros(SinularIdx.size))
        SecondFlip = False
    else:
        SinularIdx = np.where(np.diff(UpIdx) == 2)[0]
        SinularIdx = np.int32(UpIdx[SinularIdx] + np.ones(SinularIdx.size))
        np.put(Klusters, SinularIdx, np.ones(SinularIdx.size))
        SecondFlip = True
    if SecondFlip:
        SinularIdx = np.where(np.diff(DownIdx) == 2)[0]
        SinularIdx = np.int32(DownIdx[SinularIdx] + np.ones(SinularIdx.size))
        np.put(Klusters, SinularIdx, np.zeros(SinularIdx.size))
    else:
        SinularIdx = np.where(np.diff(UpIdx) == 2)[0]
        SinularIdx = np.int32(UpIdx[SinularIdx] + np.ones(SinularIdx.size))
        np.put(Klusters, SinularIdx, np.ones(SinularIdx.size))
    ###########222#######
    ###
    ###########222#######
    K1ratio = (DpwrNormal[Klusters == 1] / TpwrNormal[Klusters == 1]).mean()
    K0ratio = (DpwrNormal[Klusters == 0] / TpwrNormal[Klusters == 0]).mean()
    if K1ratio < K0ratio:
        print 'Cluster labels fliped!!!'
        Klusters = np.int8(np.logical_not(Klusters))
    k1 = Klusters.copy()
    ###
    # from now on SWS episodes are Kluster ==1
    ###
    swsIdx = np.where(Klusters == 1)[0]
    dupletIdx = np.where(np.diff(swsIdx) == 3)[0]
    dupletIdx1st = np.int32(swsIdx[dupletIdx] + np.ones(dupletIdx.size))
    dupletIdx2nd = np.int32(swsIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
    ####Repeat the same procedure for REM
    remIdx = np.where(Klusters == 0)[0]
    dupletIdx = np.where(np.diff(remIdx) == 3)[0]
    dupletIdx1st = np.int32(remIdx[dupletIdx] + np.ones(dupletIdx.size))
    dupletIdx2nd = np.int32(remIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx1st, np.zeros(dupletIdx.size))
    np.put(Klusters, dupletIdx2nd, np.zeros(dupletIdx2nd.size))
    ##########333############
    ###Now chunks of length 3!
    ##########333############
    swsIdx = np.where(Klusters == 1)[0]
    dupletIdx = np.where(np.diff(swsIdx) == 4)[0]
    dupletIdx1st = np.int32(swsIdx[dupletIdx] + np.ones(dupletIdx.size))
    dupletIdx2nd = np.int32(swsIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
    dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
    np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
    ####Repeat the same procedure for REM
    remIdx = np.where(Klusters == 0)[0]
    dupletIdx = np.where(np.diff(remIdx) == 4)[0]
    dupletIdx1st = np.int32(remIdx[dupletIdx] + np.ones(dupletIdx.size))
    try:
        dupletIdx2nd = np.int32(remIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
        dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
    except IndexError:
        print 'It\'s a bit marginal!'
    np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
    np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
    ##########444############
    ###Now chunks of length 4!
    ##########444############
    swsIdx = np.where(Klusters == 1)[0]
    dupletIdx = np.where(np.diff(swsIdx) == 5)[0]
    dupletIdx1st = np.int32(swsIdx[dupletIdx] + np.ones(dupletIdx.size))
    dupletIdx2nd = np.int32(swsIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
    dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
    dupletIdx4th = np.int32(swsIdx[dupletIdx] + 4 * np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
    np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
    np.put(Klusters, dupletIdx4th, np.ones(dupletIdx4th.size))
    ####Repeat the same procedure for REM
    remIdx = np.where(Klusters == 0)[0]
    dupletIdx = np.where(np.diff(remIdx) == 5)[0]
    dupletIdx1st = np.int32(remIdx[dupletIdx] + np.ones(dupletIdx.size))
    try:
        dupletIdx2nd = np.int32(remIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
        dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
        dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 4 * np.ones(dupletIdx.size))
    except IndexError:
        print 'It\'s a bit marginal!'
    np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
    np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
    np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
    np.put(Klusters, dupletIdx4th, np.ones(dupletIdx4th.size))
    UpIdx = np.where(Klusters == 1)[0]
    DownIdx = np.where(Klusters == 0)[0]
    SingletIdxU = np.where(np.diff(DownIdx) == 2)[0]
    SingletIdxD = np.where(np.diff(UpIdx) == 2)[0]
    while SingletIdxU.size or SingletIdxD.size:
        print SingletIdxU.size, SingletIdxD.size
        if UpIdx.size > DownIdx.size:
            SingletIdxU = np.int32(DownIdx[SingletIdxU] + np.ones(SingletIdxU.size))
            np.put(Klusters, SingletIdxU, np.zeros(SingletIdxU.size))
            SecondFlip = False
        else:
            SingletIdxD = np.where(np.diff(UpIdx) == 2)[0]
            SingletIdxD = np.int32(UpIdx[SingletIdxD] + np.ones(SingletIdxD.size))
            np.put(Klusters, SingletIdxD, np.ones(SingletIdxD.size))
            SecondFlip = True
        if SecondFlip:
            SingletIdxU = np.where(np.diff(DownIdx) == 2)[0]
            SingletIdxU = np.int32(DownIdx[SingletIdxU] + np.ones(SingletIdxU.size))
            np.put(Klusters, SingletIdxU, np.zeros(SingletIdxU.size))
        else:
            SingletIdxD = np.where(np.diff(UpIdx) == 2)[0]
            SingletIdxD = np.int32(UpIdx[SingletIdxD] + np.ones(SingletIdxD.size))
            np.put(Klusters, SingletIdxD, np.ones(SingletIdxD.size))
        del SecondFlip
        ###########222#######
        ###
        ###########222#######
        K1ratio = (DpwrNormal[Klusters == 1] / TpwrNormal[Klusters == 1]).mean()
        K0ratio = (DpwrNormal[Klusters == 0] / TpwrNormal[Klusters == 0]).mean()
        if K1ratio < K0ratio:
            print 'Cluster labels fliped!!!'
            Klusters = np.int8(np.logical_not(Klusters))
        k1 = Klusters.copy()
        ###
        # from now on SWS episodes are Kluster ==1
        ###
        swsIdx = np.where(Klusters == 1)[0]
        dupletIdx = np.where(np.diff(swsIdx) == 3)[0]
        dupletIdx1st = np.int32(swsIdx[dupletIdx] + np.ones(dupletIdx.size))
        dupletIdx2nd = np.int32(swsIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
        ####Repeat the same procedure for REM
        remIdx = np.where(Klusters == 0)[0]
        dupletIdx = np.where(np.diff(remIdx) == 3)[0]
        dupletIdx1st = np.int32(remIdx[dupletIdx] + np.ones(dupletIdx.size))
        dupletIdx2nd = np.int32(remIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx1st, np.zeros(dupletIdx.size))
        np.put(Klusters, dupletIdx2nd, np.zeros(dupletIdx2nd.size))
        ##########333############
        ###Now chunks of length 3!
        ##########333############
        swsIdx = np.where(Klusters == 1)[0]
        dupletIdx = np.where(np.diff(swsIdx) == 4)[0]
        dupletIdx1st = np.int32(swsIdx[dupletIdx] + np.ones(dupletIdx.size))
        dupletIdx2nd = np.int32(swsIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
        dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
        np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
        ####Repeat the same procedure for REM
        remIdx = np.where(Klusters == 0)[0]
        dupletIdx = np.where(np.diff(remIdx) == 4)[0]
        dupletIdx1st = np.int32(remIdx[dupletIdx] + np.ones(dupletIdx.size))
        try:
            dupletIdx2nd = np.int32(remIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
            dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
        except IndexError:
            print 'It\'s a bit marginal!'
        np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
        np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
        ##########444############
        ###Now chunks of length 4!
        ##########444############
        swsIdx = np.where(Klusters == 1)[0]
        dupletIdx = np.where(np.diff(swsIdx) == 5)[0]
        dupletIdx1st = np.int32(swsIdx[dupletIdx] + np.ones(dupletIdx.size))
        dupletIdx2nd = np.int32(swsIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
        dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
        dupletIdx4th = np.int32(swsIdx[dupletIdx] + 4 * np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
        np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
        np.put(Klusters, dupletIdx4th, np.ones(dupletIdx4th.size))
        ####Repeat the same procedure for REM
        remIdx = np.where(Klusters == 0)[0]
        dupletIdx = np.where(np.diff(remIdx) == 5)[0]
        dupletIdx1st = np.int32(remIdx[dupletIdx] + np.ones(dupletIdx.size))
        try:
            dupletIdx2nd = np.int32(remIdx[dupletIdx] + 2 * np.ones(dupletIdx.size))
            dupletIdx3rd = np.int32(swsIdx[dupletIdx] + 3 * np.ones(dupletIdx.size))
            dupletIdx4th = np.int32(swsIdx[dupletIdx] + 4 * np.ones(dupletIdx.size))
        except IndexError:
            print 'It\'s a bit marginal!'
        np.put(Klusters, dupletIdx1st, np.ones(dupletIdx.size))
        np.put(Klusters, dupletIdx2nd, np.ones(dupletIdx2nd.size))
        np.put(Klusters, dupletIdx3rd, np.ones(dupletIdx3rd.size))
        np.put(Klusters, dupletIdx4th, np.ones(dupletIdx4th.size))

        UpIdx = np.where(Klusters == 1)[0]
        DownIdx = np.where(Klusters == 0)[0]
        SingletIdxU = np.where(np.diff(DownIdx) == 2)[0]
        SingletIdxD = np.where(np.diff(UpIdx) == 2)[0]
        ##############################################################
        #####Creating sleep state(REM/SWS) dataframe
        swsInit = np.where(np.diff(Klusters) > 0)[0]
        remInit = np.where(np.diff(Klusters) < 0)[0]
        if swsInit.size != remInit.size:
            if not swsInit[0]:
                swsInit = swsInit[1:]
            if not remInit[0]:
                remInit = remInit[1:]
            if swsInit[-1] == Klusters.size - 2:
                swsInit = swsInit[:-1]
            if remInit[-1] == Klusters.size - 2:
                remInit = remInit[:-1]
        print swsInit.size, remInit.size
        swsPeriods = []
        remPeriods = []
        if remInit[0] < swsInit[0]:
            swsPeriods.append([0, remInit[0]])
            for ii in range(remInit.size):
                try:
                    remPeriods.append([remInit[ii], swsInit[ii]])
                    swsPeriods.append([swsInit[ii], remInit[ii + 1]])
                except IndexError:
                    a = 1
        else:
            remPeriods.append([0, swsInit[0]])
            for ii in range(remInit.size):
                try:
                    swsPeriods.append([swsInit[ii], remInit[ii]])
                    remPeriods.append([remInit[ii], swsInit[ii + 1]])
                except IndexError:
                    a = 1
        swsPeriods = np.array(swsPeriods)
        remPeriods = np.array(remPeriods)
        ######Saving Immobility dataFrame
        SleepTimeIndex.to_pickle(animalPath + animalID + '/' +
                                 '_' + sessionDate + '_' + sessionID + '.df')
        ksum = 0
        pathSplit = lfp.tags['file'].split('/')
        animalcnd = np.array(['Rat' in item for item in pathSplit])
        datecnd = np.array(['LinearTrack' in item for item in pathSplit])
        sessioncnd = np.array(['-sleep' in item for item in pathSplit])
        animal = np.array(pathSplit)[animalcnd][0]
        date = np.array(pathSplit)[datecnd][0].split('_')[1]
        session = np.array(pathSplit)[sessioncnd][0].split('-')[-1]
        print session, date, animal
        ###First Sleep period and initiating Dataframe
        validSWSepisodes = np.where(np.logical_and(swsPeriods[:, 0] > ksum,
                                                   swsPeriods[:, 1] < ksum +
                                                   SleepTimeIndex.durationSeconds[0]))[0]
        validREMepisodes = np.where(np.logical_and(remPeriods[:, 0] > ksum,
                                                   remPeriods[:, 1] < ksum +
                                                   SleepTimeIndex.durationSeconds[0]))[0]
        SleepPeriodsDF = pd.DataFrame([[SleepTimeIndex.initTime[0] + 1e3 * (swsPeriods[0, 0] - ksum),
                                        SleepTimeIndex.initTime[0] + 1e3 * (swsPeriods[0, 1] - ksum),
                                        'sws', session, date, animal]],
                                      columns=['t0', 't1', 'epoch', 'session', 'date', 'animal'])
        for item in validSWSepisodes[1:]:
            df = pd.DataFrame([[SleepTimeIndex.initTime[0] + 1e3 * (swsPeriods[item, 0] - ksum),
                                SleepTimeIndex.initTime[0] + 1e3 * (swsPeriods[item, 1] - ksum),
                                'sws', session, date, animal]],
                              columns=['t0', 't1', 'epoch', 'session', 'date', 'animal'])
            SleepPeriodsDF = SleepPeriodsDF.append(df, ignore_index=True)
            del df
        print validSWSepisodes
        print validREMepisodes
        for item in validREMepisodes:
            df = pd.DataFrame([[SleepTimeIndex.initTime[0] + 1e3 * (remPeriods[item, 0] - ksum),
                                SleepTimeIndex.initTime[0] + 1e3 * (remPeriods[item, 1] - ksum),
                                'rem', session, date, animal]],
                              columns=['t0', 't1', 'epoch', 'session', 'date', 'animal'])
            SleepPeriodsDF = SleepPeriodsDF.append(df, ignore_index=True)
            del df
        ksum += SleepTimeIndex.durationSeconds[0]
        for ii in SleepTimeIndex.index[1:]:
            # fig,ax = pl.subplots(1,1)
            validSWSepisodes = np.where(np.logical_and(swsPeriods[:, 0] > ksum,
                                                       swsPeriods[:, 1] < ksum +
                                                       SleepTimeIndex.durationSeconds[ii]))[0]
            validREMepisodes = np.where(np.logical_and(remPeriods[:, 0] > ksum,
                                                       remPeriods[:, 1] < ksum +
                                                       SleepTimeIndex.durationSeconds[ii]))[0]
            for item in validSWSepisodes[1:]:
                df = pd.DataFrame([[SleepTimeIndex.initTime[0] + 1e3 * (swsPeriods[item, 0] - ksum),
                                    SleepTimeIndex.initTime[0] + 1e3 * (swsPeriods[item, 1] - ksum),
                                    'sws', session, date, animal]],
                                  columns=['t0', 't1', 'epoch', 'session', 'date', 'animal'])
                SleepPeriodsDF = SleepPeriodsDF.append(df, ignore_index=True)
                del df
            for item in validREMepisodes:
                df = pd.DataFrame([[SleepTimeIndex.initTime[0] + 1e3 * (remPeriods[item, 0] - ksum),
                                    SleepTimeIndex.initTime[0] + 1e3 * (remPeriods[item, 1] - ksum),
                                    'rem', session, date, animal]],
                                  columns=['t0', 't1', 'epoch', 'session', 'date', 'animal'])
                SleepPeriodsDF = SleepPeriodsDF.append(df, ignore_index=True)
                del df
            # for item in validREMepisodes:
            #    ax.axvspan(remPeriods[item,0]-ksum,remPeriods[item,1]-ksum,alpha=0.13,color='c')

            # print validSWSepisodes
            ksum += SleepTimeIndex.durationSeconds[ii]
cd
            ######PLOT CHECK
            fig, ax = pl.subplots(1, 1, figsize=[30, 10])
            ax.plot(lfp.timeAxis, zScore(lfp.signal), lw=1, alpha=0.82)
            ax.plot(1e3 * ts, zScore(r), alpha=0.92, lw=1)
            # for item in sleepList:
            #    ax.axvspan(1e3*item[0],1e3*item[1],alpha=0.5,color='r')
            for index, row in SleepTimeIndex.iterrows():
                ax.axvspan(row.initTime, row.initTime + row.duration, alpha=0.5, color='c')
            for idx in SleepPeriodsDF.index:
                item = SleepPeriodsDF.loc[idx]
                if (item.epoch == 'sws'):
                    ax.axvspan(item.t0, item.t1, alpha=0.25, color='b')
                else:
                    ax.axvspan(item.t0, item.t1, alpha=0.25, color='r')