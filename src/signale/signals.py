"""
signale.signals
===============

A module for signal classes.
"""
__author__ = ("KT")
__version__ = "4.1, August 2013"


# python modules
import inspect


# other modules
import numpy as np

from scipy.optimize import fminbound
from scipy.special import erf
from scipy.stats import norm
from scipy.signal import butter,filtfilt, freqz
import matplotlib.pyplot as pl
import matplotlib as mpl


# custom made modules
import custom_plot

# package modules
# #import io



###################################################### FUNCTIONS



def ExpKer(length,tau):
    '''A function to get a causal exponential kernel.
    
    math: k(t) = H(t) \exp{\left(- \frac{t}{\tau_R}\right)}
          H(t) = 0 for x< 0 and 1 otherwise
    '''

    
    x = np.linspace(-1,1,length)
    x = np.piecewise(x,[x < 0, x>=0],[0,1])
    kernel = x * np.exp(-1 *np.linspace(-1,1,length) / tau)
    return kernel

def fitGauss(x, data):
    """ Fit a gaussian to data.

    Using method of moments.

    Paramters:
    x ... x axis of data
    data ... 1D array of data

    Returns:
    gauss ... fitted gaussian
    """

    # X = np.arange(data.size)
    x_bar = np.sum(x * data) / np.sum(data)
    width = np.sqrt(np.abs(np.sum((x - x_bar) ** 2 * data) / np.sum(data)))
    amp = data.max()

    gauss = norm.pdf(x, x_bar, width)
    gauss /= gauss.max()
    gauss *= amp

    return gauss

def butter_bandpass(lowcut, highcut, fs, order=5):
    """
    A function to design a Butterworth bandpass filter.
    Parameters
        ----------
        sig : This is the first signal.
        lowcut:
        highcut:
        fs: Sampling frequency.
                
        Returns
        ----------
        b:
        a:
        
        See also
        ----------
        scipy.signal
        
        Notes
        ----------
    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(sig, lowcut, highcut, fs, order=5 ):
    '''
    A function to perform a bandpass filter on "data".
    Parameters
        ----------
        sig : This is the first signal.
        lowcut:
        highcut:
        fs: Sampling frequency.
        order:
        Returns
        ----------
        filteredSig
        
        See also
        ----------
        scipy.signal
        
        Notes
        ----------
        
    '''
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    filteredSig = filtfilt(b, a, sig)
    return filteredSig

def freq_response(lowcut,highcut,fs,order= [5],fig=None, ax=None,display=None):
    """
    A function to plot the response function of the filter.
    
    Parameters
    ----------
    lowcut:
    highcut:
    fs: Sampling frequency.
    order:
    
    Returns
    ----------
    responseCurve
           
    See also
    ----------
    scipy.signal
        
    Notes
    ----------
    """
    if not fig:
        fig = pl.figure(figsize=(12, 7))
    if not ax:
        ax = fig.add_subplot(111)
    for item in order:
           b, a = butter_bandpass(lowcut, highcut, fs, order=item)
           w, h = freqz(b, a, worN=2000)
           pl.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % item)
   
    pl.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],'--',label="sqrt(0.5)")
    pl.xlabel('Frequency (Hz)')
    pl.ylabel('Gain')
    pl.grid(True)
    pl.legend(loc='best')
    pl.xlim(0.75*lowcut,1.5*highcut) 
    if display:
        pl.show()
    return fig,ax



def pairs(a):
    """ Generator function yielding all item pairs of a list.

    Parameters:
    a ... list to iterate over

    """

    l = len(a)
    for i in range(l):
        for j in range(i + 1, l):
            yield a[i], a[j]


def rcc(lindat, phidat, display=False):
    """
    Calculate slope and correlation value on circular data (phidat).

    phidat is expected in radiants.
    """
    global aopts, Rs
    Rs = []
    aopts = []

    phidat = np.float_(phidat)  # make lindat to floats
    lindat = np.float_(lindat)  # make lindat to floats
   # lindat = lindat/lindat.max()    # normalize lindat


    # bounds of slope
    abound = [-2., 2.]

    # starting points of maximization
    Nrep = 20
    da = abs(abound[1] - abound[0]) / Nrep
    astart = min(abound)

    Rfunc = lambda a:-np.absolute(np.mean(np.exp(1j * (phidat - 2.*np.pi * a * lindat))))

    aopt = np.nan
    R = -10 ** 10
    for n in range(Nrep):
        a0 = astart + n * da;
        a1 = astart + (n + 1) * da;

        returnValues = fminbound(Rfunc, a0, a1, full_output=True, disp=1)  # xtol=1e-10
        # print returnValues
        # atmp = returnValues[0][0]
        atmp = returnValues[0]
        rtmp = returnValues[1]

        if display:
            aopts.append(atmp)
            Rs.append(-rtmp)

        if -rtmp > R:
         #   print rtmp, R, aopt, atmp
            R = -rtmp
            aopt = atmp

    # phase offset
    v = np.mean(np.exp(1j * (phidat - 2 * np.pi * aopt * lindat)))
    phi0 = np.angle(v)

    theta = np.angle(np.exp(2 * np.pi * 1j * np.abs(aopt) * lindat))

    thmean = np.angle(np.sum(np.exp(1j * theta)))
    phmean = np.angle(np.sum(np.exp(1j * phidat)))

    sthe = np.sin(theta - thmean)
    sphi = np.sin(phidat - phmean)

    c12 = np.sum(sthe * sphi)
    c11 = np.sum(sthe * sthe)
    c22 = np.sum(sphi * sphi)

    rho = c12 / np.sqrt(c11 * c22)

    lam22 = np.mean(sthe ** 2.*sphi ** 2)
    lam20 = np.mean(sphi ** 2)
    lam02 = np.mean(sthe ** 2)
    tval = rho * np.sqrt(lindat.size * lam20 * lam02 / lam22)

    p = 1 - erf(tval / np.sqrt(2))

    if display:
        fig = pl.figure(1, figsize=(10, 10))
        ax = fig.add_subplot(2, 2, 1)
        ax.plot(lindat, phidat, '.k')
        ax.plot(lindat, phidat - 2 * np.pi, '.k')
        ax.plot(lindat, phi0 + 2 * np.pi * aopt * lindat, 'r-')

        ax = fig.add_subplot(2, 2, 2)
        ax.plot(lindat, theta, '.k')

        ax = fig.add_subplot(2, 3, 4)
        ax.plot(sthe, sphi, '.b')

        ax = fig.add_subplot(2, 3, 5)
        ax.hist(sthe * sphi, 20)

        ax = fig.add_subplot(2, 3, 6)
        ax.plot(aopts, Rs, '.-')
        ax.set_xlabel('a')
        ax.set_ylabel('R')

        pl.show()

    return rho, p, R, aopt, phi0
def fft(sig,sampFreq=None):
    '''
    Calculates the FFT of a signal(sig) with sampling frequency sampFreq.
       
        Parameters
        ----------
        sig : This is the first signal.
        sampFreq: Sampling frequency.
                
        Returns
        ----------
        freq:
        spPower:
        
        See also
        ----------
        signale.cscs.NeuralynxCSC.ripple_recorder, np.where
        
        Notes
        ----------
        

     '''
    sig = np.array(sig)
    nsig = sig.size
    sp = np.fft.rfft(sig)
    spPower = 2 * sp * np.conj(sp)
    freq = np.fft.fftfreq(nsig, d = 1./sampFreq)[0:nsig/2+1]
    return freq,spPower

def filter(sig,freq=np.array([]),sp = np.array([]),minFreq = None,maxFreq=None,sampFreq=None):
        '''
            Calculates the bandPassed filtered of a signal(sig) between min_freq and max_freq frequencies.
       
        Parameters
        ----------
        sig : This is the first signal.
        sampFreq: Sampling frequency.
                
        Returns
        ----------
        signal_filtered: 
        
        See also
        ----------
        signale.cscs.NeuralynxCSC.ripple_recorder, np.where
        
        Notes
        ----------
        

         '''
        if not sp.any():
            if not sampFreq:
                print "Please enter the sampling frequency!!!"
            else:
                freq,sp = fft(sig,sampFreq)
        if not minFreq:
            print "Please enter the minimum frequency!"
        if not maxFreq:
                print "Please enter the maximum frequency!"

                # band pass filtering, zero phase
        sp_filtered = np.zeros_like(sp)
        for i, f in enumerate(freq):
            if abs(f) >= minFreq and abs(f) <= maxFreq:
                sp_filtered[i] = sp[i]

        # backtransform filtered signal to time space
        sig_filtered = np.fft.irfft(sp_filtered)
        sig_filtered *= sig.size  # rescale from normalization of the fourier components
        sig_filtered = sig_filtered[:sig.size]
        return sig_filtered
def _getUnitFactor(timeUnit, newUnit):
    """
    Returns factor to change time unit.
    """

    factor = 1.
    if timeUnit == 'us':
        if newUnit == 'ms':
            factor /= 1000
        elif newUnit == 's':
            factor /= 1000 * 1000
        elif newUnit == 'min':
            factor /= 1000 * 1000 * 60
        elif newUnit == 'h':
            factor /= 1000 * 1000 * 60 * 60
    elif timeUnit == 'ms':
        if newUnit == 'us':
            factor *= 1000
        elif newUnit == 's':
            factor /= 1000
        elif newUnit == 'min':
            factor /= 1000 * 60
        elif newUnit == 'h':
            factor /= 1000 * 60 * 60
    elif timeUnit == 's':
        if newUnit == 'us':
            factor *= 1000 * 1000
        elif newUnit == 'ms':
            factor *= 1000
        elif newUnit == 'min':
            factor /= 60
        elif newUnit == 'h':
            factor /= 60 * 60
    elif timeUnit == 'min':
        if newUnit == 'us':
            factor *= 60 * 1000 * 1000
        elif newUnit == 'ms':
            factor *= 60 * 1000
        elif newUnit == 's':
            factor *= 60
        elif newUnit == 'h':
            factor /= 60
    elif timeUnit == 'h':
        if newUnit == 'us':
            factor *= 60 * 60 * 1000 * 1000
        elif newUnit == 'ms':
            factor *= 60 * 60 * 1000
        elif newUnit == 's':
            factor *= 60 * 60
        elif newUnit == 'min':
            factor *= 60

    return factor





###################################################### CLASSES

class signal(object):

    def __init__(self):
        self.meta = {}
        self.tags = {}  # for comments/tags to the spiketrains
        self.timeUnit = 'ms'

    def addTags(self, id, **kwargs):
        """
        For adding and changing comments/tags to the signal object.
        """
        pass

    def showTags(self):
        """
        For print the comments/tags to the signal object on the command line.
        """
        for key in self.tags:
            print key, ':', self.tags[key]
        print ''

    def changeTimeUnit(self):
        print 'WARNING: Not implemented!'




class NeuralynxEvents(signal):
    """
    A class for Neuralynx event data.
    """

    def __init__(self, times, eventStrings):
        signal.__init__(self)
        self.times = np.array(times)
        self.eventStrings = eventStrings

        self.t_start = -1.
        self.t_stop = -1.
        self.t_startRecording = -1.
        self.t_stopRecording = -1.

        self.t_startRecording = self.times[0]  # remember the time of the start of the recording separatly
        if self.eventStrings[0] == 'Starting Recording':
            self.times = self.times[1:]  # remove entry of the start of the recording
            self.eventStrings = self.eventStrings[1:]  # remove entry of the start of the recording
        # if self.eventStrings[0].find('Digital Lynx Parallel Input Port TTL')+1:
        self.t_start = self.times[0]  # remember the time of the first TTL pulse as t_start

        self.t_stopRecording = self.times[-1]  # remember the time of the end of the recording separatly
        if self.eventStrings[-1] == 'Stopping Recording':
            self.times = self.times[:-1]  # remove corresponding entry
            self.eventStrings = self.eventStrings[:-1]  # remove corresponding entry

    def changeTimeUnit(self, newUnit='ms'):

        factor = _getUnitFactor(self.timeUnit, newUnit)

        # change the times
        self.times *= factor
        self.t_start *= factor
        self.t_stop *= factor
        self.timeUnit = newUnit

        self.t_startRecording *= factor
        self.t_stopRecording *= factor


    def purge(self, numItems):
        """
        Cut away numItems.
        """
        self.times = self.times[numItems:]
        self.eventStrings = self.eventStrings[numItems:]


    def time_offset(self, offset):
        """
        Add a time offset to the NeuralynxEvents object. times, t_start and t_stop are
        shifted from offset.
        """

        if not self.times.size:
            print 'WARNING: times array is empty, can not offset time.'
            return

        self.times += offset
        self.t_start = self.times[0]
        self.t_stop = self.times[-1]

        self.t_startRecording += offset
        self.t_stopRecording += offset


