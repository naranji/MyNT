'''
Created on 26 Mar 2014
This module contains the methods for comparing and distinguishing the REM, SWS periods during sleep.

@author: chenani
'''
from __future__ import division
import os, sys
import numpy as np
import matplotlib as mpl
mpl.matplotlib_fname()
import matplotlib.pyplot as pl
import pickle as pkl
import scipy.signal as scsig
from matplotlib.backends.backend_pdf import PdfPages
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
# if len(sys.argv) > 1:
#     try: 
#         fileName = sys.argv[1]
#     except:
#         print "Usage:", sys.argv[0], "foo.ephys"; sys.exit(1)
# else:
#     fileName = raw_input('Enter command line arguments: ').split()[0]
# rec = pkl.load(open(fileName, "rb"))
# csc = rec.csc
# FUNCTIONS
def rem_content(x =signale.cscs.NeuralynxCSC, threshold = None):
    if not hasattr(x,'th_del_ratio'):
        x.theta_delta_ratio()
    if threshold:
        x.REM_detector(thresh = threshold)
    else:
        x.REM_detector(thresh = x.th_del_ratio[1].mean()+x.th_del_ratio[1].std())
    x.SWS_signal()
    x.REM_signal()
    sws = []
    rem = []
    for item in x.sws_signal:
        for jtem in item:
            sws.append(jtem)
    for item in x.rem_signal:
        for jtem in item:
            rem.append(jtem)
    sws = np.array(sws)
    rem = np.array(rem)
    return sws,rem 
def steppify(arr,isX=False,interval=0):
    """
    Converts an array to double-length for step plotting
    """
    if isX and interval==0:
        interval = abs(arr[1]-arr[0]) / 2.0
        newarr = np.array(zip(arr-interval,arr+interval)).ravel()
        return newarr
def downsampler(arr,factor):
    newArr = np.array([])
    for item in arr:
        if np.mod(np.where(arr==item)[0],factor) == 0:
            newArr = np.append(newArr,item)
    return newArr
def fft(x,y,sampleFreq = None, display = False,fig=None,ax=None,step=False):
    '''
    This function makes and compares the power spectrum of two signals. This is mainly designed for comparison between SWS and REM
    signals extracted from sleep LFP signal.
       
        Parameters
        ----------
        x : This is the first signal.
        y : the second signal.
        sampleFreq: Sampling frequency.
        display: boolean variable, you can switch on if you like to see the plots.
        
        Returns
        ----------
        
        
        See also
        ----------
        signale.cscs.NeuralynxCSC.ripple_recorder, numpy.where
        
        Notes
        ----------
        

     ''' 
    x = np.array(x)
    y = np.array(y)
    nx = x.size
    ny = y.size
    xsp = np.fft.rfft(x)
    ysp = np.fft.rfft(y)
    xspPower = 2 * xsp * np.conj(xsp)
    yspPower = 2 * ysp * np.conj(ysp)
    xfreq = np.fft.fftfreq(nx, d = 1./sampleFreq)[0:nx/2+1]
    yfreq = np.fft.fftfreq(ny, d = 1./sampleFreq)[0:ny/2+1]
    xwin = scsig.gaussian(560, 80)
    ywin = scsig.gaussian(560,80)
    xtc = np.where(xfreq > 80)[0][0]
    ytc = np.where(yfreq > 80)[0][0]
    hehe = np.convolve(xspPower, xwin,"same" )
    hoho = np.convolve(yspPower, ywin,"same")
    hehe /= hehe[xtc:].max()
    hoho /= hoho[ytc:].max()
    idx = (np.abs(yfreq-150)).argmin()
    if display:
        with PdfPages('/home/chenani/fft_plots.pdf') as pdf:
            if not fig:
                fig = pl.figure()
            if not ax:
                ax = fig.add_subplot(111)
            if step:
                xfdown = downsampler(xfreq,80)
                hedown = downsampler(hehe,80)
                yfdown = downsampler(yfreq,80)
                hodown = downsampler(hoho,80)
                ax.bar(xfdown,hedown,width=xfdown[1]-xfdown[0],color='green',linewidth=1,label = 'real')
                ax.bar(yfdown,hodown,width=yfdown[1]-yfdown[0],color='blue',linewidth=1,label='VR')
                #ax.step(xfdown,hedown)
                #ax.fill_between(steppify(xfdown[1:xfdown.size],isX=True),steppify(hedown[1:hedown.size],isX=True) * 0,steppify(hedown[1:hedown.size],isX=True),facecolor='green',alpha=0.6)
            else:
                ax.fill_between(xfreq,hehe,facecolor='green',alpha=0.6)
                ax.fill_between(yfreq,hoho,facecolor='blue',alpha=0.6)
            ax.set_ylabel("Power(normalized)")
            pl.xlabel("Frequency(Hz)")
            pl.xlim(120,220)
            pl.ylim(0,2*hoho[idx])
            pl.legend(loc='upper right')
#         ax[1].plot(xfreq,xspPower,'r',alpha = 0.6)
#         ax[1].plot(yfreq,yspPower,'r',alpha = 0.4)
#         ax[1].set_ylabel("Power")
            
        
            pdf.savefig()
        return xfdown,hedown,yfdown,hodown
    else : 
        return hehe,xfreq,hoho,yfreq
        
        