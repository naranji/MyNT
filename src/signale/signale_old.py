"""
A module for signal classes. 
"""
__author__ = ("KT", "Moritz Dittmeyer", "Franziska Hellmundt", "Christian Leibold")
__version__ = "3.2, February 2013"


# python modules
import inspect, struct, os, sys


# other modules
import numpy

import scipy, scipy.signal
from scipy.optimize import fminbound 
from scipy.special import erf
from scipy.stats import norm

import NeuroTools.signals as NTsig

import matplotlib.pyplot as pl
import matplotlib.mlab
import matplotlib as mpl
import matplotlib.nxutils as mpl_nxutils


# custom made modules
import custom_plot

###################################################### FUNCTIONS



def findMaximaMinima(a, findMaxima=1, findMinima=1, or_equal=1):
    """ Find the maxima and minima in an 1D array.
    
    Paramters:
    a ... 1D array
    findMaxima, findMinima ... boolean, if ought to look for maxima and/or minima [defaults 1]
    
    Returns:
    dictionary{'maxima_number', 'minima_number', 'maxima_indices', 'minima_indices'}
    """
    gradients = numpy.diff(a)

    maxima_num = 0
    minima_num = 0
    max_indices = []
    min_indices = []    

    if or_equal:  # includes only the next entries and they can also be equal to zero
        for i, g in enumerate(gradients[:-1]):
            if findMaxima and cmp(g, 0) >= 0 and cmp(gradients[i + 1], 0) <= 0 and g != gradients[i + 1]:
                maxima_num += 1
                max_indices.append(i + 1)     

            if findMinima and cmp(g, 0) <= 0 and cmp(gradients[i + 1], 0) >= 0 and g != gradients[i + 1]:
                minima_num += 1
                min_indices.append(i + 1)
    else:  # includes also the previous and next entries and also expects them not to be equal to zero
        for i, g in enumerate(gradients[1:-1]):
            j = i + 1  # get real index
            if findMaxima and cmp(g, 0) >= 0 and cmp(gradients[j - 1], 0) > 0 and cmp(gradients[j + 1], 0) < 0:
                maxima_num += 1
                max_indices.append(j + 1)     

            if findMinima and cmp(g, 0) <= 0 and cmp(gradients[j - 1], 0) < 0 and cmp(gradients[j + 1], 0) > 0:
                minima_num += 1
                min_indices.append(j + 1)
   
    return {'maxima_number': maxima_num , 'minima_number': minima_num, \
            'maxima_indices': numpy.array(max_indices), 'minima_indices': numpy.array(min_indices)}  


def findNearest(array, value):
    """ Find the entry nearest to a value in an 1D array.
    
    Paramters:
    array ... 1D array
    value ... value to find nearest entry for
    
    Returns:
    a tuple (index, array[index])
    """
    index = (numpy.abs(array - value)).argmin()
    return (index, array[index])



def fitGauss(x, data):
    """ Fit a gaussian to data.
    
    Using method of moments.
    
    Paramters:
    x ... x axis of data
    data ... 1D array of data
    
    Returns:
    gauss ... fitted gaussian
    """
    
    # X = numpy.arange(data.size)
    x_bar = numpy.sum(x * data) / numpy.sum(data)
    width = numpy.sqrt(numpy.abs(numpy.sum((x - x_bar) ** 2 * data) / numpy.sum(data)))
    amp = data.max()

    gauss = norm.pdf(x, x_bar, width)
    gauss /= gauss.max()
    gauss *= amp
    
    return gauss


def load_tFile(fileName, showHeader=False):
    """
    For reading AD Redishs's MClust 3.5 t-files with python/numpy.
    These files contain events/spikes data.
    
    Returns a Neurotools.signals SpikeTrain object. 
    """
    spikes = []
    header = ''
    with open(fileName, 'rb') as f:
        headerFinished = False
        while not headerFinished:
            line = f.readline()
            if line.startswith('%%ENDHEADER'):
                headerFinished = True
            header += line
        while True:
            chunk = f.read(4)
            if chunk:
                spike = struct.unpack('>L', chunk)[0]  # data in MClust t files is stored as UInt32 big endian
                spikes.append(spike)
            else:
                break
    spikes = placeCell_spikezug(numpy.float_(spikes) / 10.)  # convert to ms, 
                                                               # time in MClust t-files is in tens of ms
    
    if showHeader:
        print header

    return spikes


def load_ncsFile(fileName, showHeader=False):
    """
    For loading Neuralynx continuously sampled channel (CSC) recorded data
    with python/numpy.
    Returns a NeuralynxCSC object. 
    """
    timeStamps = []
    dataSamples = []
    dataSamplesNums = []
    
    with open(fileName, 'rb') as f:
    
        header = f.read(16 * 2 ** 10)  # the first 16 kB are an ASCII text header
        
        if showHeader:
            print header

        count = -1
        while True:
            qwTimeStamp = f.read(8)  # UInt64, Cheetah timestamp for this record. This value is in microseconds.
            dwChannelNumber = f.read(4)  # UInt32, channel number for this record.
            dwSampleFreq = f.read(4)  # UInt32, sampling frequency reported by the acquisition hardware when recording this record.
            dwNumValidSamples = f.read(4)  # UInt32, Number of values in snSamples containing valid data.
            snSamples = f.read(512 * 2)  # Int16[], Data points for this record.

            if qwTimeStamp:
                qwTimeStamp = struct.unpack('L', qwTimeStamp)[0]
                dwChannelNumber = struct.unpack('I', dwChannelNumber)[0]
                dwSampleFreq = struct.unpack('I', dwSampleFreq)[0]
                dwNumValidSamples = struct.unpack('I', dwNumValidSamples)[0]
                
                channel = dwChannelNumber
                sampleFreq = dwSampleFreq
                validSamples = dwNumValidSamples
                
                snSamples = [snSamples[i:i + 2] for i in range(0, snSamples.__len__(), 2)]
                samples = []
                for sample in snSamples:
                    sample = struct.unpack('h', sample)[0]
                    samples.append(sample)
                snSamples = numpy.array(samples)
                
                timeStamps.append(qwTimeStamp)
                dataSamples.append(snSamples)
                dataSamplesNums.append(snSamples.size)
                
                if dataSamplesNums[-1] != validSamples:
                    print 'WARNING: Number of samples and number of valid samples not the same.'
                
                
                # print some data?
                if count < 0 and showHeader:
                    count += 1
                    print '----> count', count      
                    print qwTimeStamp         
                    print channel          
                    print sampleFreq             
                    print validSamples                   
                    # print snSamples
                    print snSamples.shape              
                    print ''
            else:
                break
    
    timeStamps = numpy.array(timeStamps) / 1000.  # change to ms, Neuralynx time stamps are in microseconds
    dataSamples = numpy.array(dataSamples)
    dataSamples = dataSamples.flatten()
    
    # integrity check, after extracellpy by Santiago Jaramillo
    if numpy.any(numpy.diff(numpy.diff(timeStamps))):
        print('WARNING: Not all records are contiguous. Packets lost?')
    
    csc = NeuralynxCSC(timeStamps, dataSamples, sampleFreq)
    csc.tags['file'] = fileName
    csc.tags['path'] = os.getcwd()
    csc.tags['channel'] = channel
    csc.tags['validSamples'] = validSamples
    csc.header = header
    csc.dataSamplesNums = dataSamplesNums
    
    csc.recalc_timeAxis()  # recalculate time values and timeAxis using the validSamples
    
    return csc


def load_nevFile(fileName, showHeader=False):
    """
    For loading Neuralynx Event files with python/numpy.
    Returns a NeuralynxEvent object. 
    """
    timeStamps = []
    eventStrings = []
    
    with open(fileName, 'rb') as f:
    
        header = f.read(16 * 2 ** 10)  # the first 16 kB are an ASCII text header
        
        if showHeader:
            print header
        
        count = -3
        while True:
            nstx = f.read(2)  # Int16, Reserved
            npkt_id = f.read(2)  # Int16, ID for the originating system of this packet.
            npkt_data_size = f.read(2)  # Int16, This value should always be two (2).
            qwTimeStamp = f.read(8)  # UInt64, Cheetah timestamp for this record. This value is in microseconds.
            nevent_id = f.read(2)  # Int16, ID value for this event.
            nttl = f.read(2)  # Int16, Decimal TTL value read from the TTL input port.
            ncrc = f.read(2)  # Int16, Record CRC check from Cheetah. Not used in consumer applications.
            ndummy1 = f.read(2)  # Int16, Reserved
            ndummy2 = f.read(2)  # Int16, Reserved
            dnExtra = f.read(8 * 4)  # Int32[], Extra bit values for this event. 
                                                # This array has a fixed length of eight (8).
            eventString = f.read(128)  # Event string associated with this event record. This string has a maximum length of 128 characters.

            if nstx:
                nstx = struct.unpack('h', nstx)[0]
                npkt_id = struct.unpack('h', npkt_id)[0]
                npkt_data_size = struct.unpack('h', npkt_data_size)[0]
                qwTimeStamp = struct.unpack('L', qwTimeStamp)[0]
                nevent_id = struct.unpack('h', nevent_id)[0]
                nttl = struct.unpack('h', nttl)[0]
                ncrc = struct.unpack('h', ncrc)[0]
                ndummy1 = struct.unpack('h', ndummy1)[0]
                ndummy2 = struct.unpack('h', ndummy2)[0]
                
                timeStamps.append(qwTimeStamp)
                eventStrings.append(eventString.split('\x00')[0])
                
                # print some data?
                if count < 0:
                    count += 1
                    print '----> count', count 
                    print nstx                    
                    print npkt_id                
                    print npkt_data_size        
                    print qwTimeStamp         
                    print nevent_id          
                    print nttl             
                    print ncrc                   
                    print ndummy1              
                    print ndummy2                
                    print dnExtra              
                    print eventString
                    print ''
            else:
                break
        
        timeStamps = numpy.array(timeStamps) / 1000.  # change to ms, Neuralynx time stamps are in microseconds
        
        nev = NeuralynxEvents(timeStamps, eventStrings)
        nev.tags['file'] = fileName
        nev.tags['path'] = os.getcwd()
        
        return nev



def load_rawFile(fileName, exclude=[], showHeader=False):
    """
    For loading RAW data file exported from Multi Channel Systems MCRack with python/numpy.
    """
    data = []
    dataSamplesNums = []
    
    header = ''
    with open(fileName, 'rb') as f:
        headerFinished = False
        while not headerFinished:
            line = f.readline()
            if line.startswith('EOH'):
                headerFinished = True
            header += line
            
            # get metadata from header
            if line.startswith('Sample rate'):
                sampleFreq = float(line.split('=')[1])
            elif line.startswith('Streams'):
                channels = line.split('=')[1]
                channels = channels.split(';')
                channels[-1] = channels[-1][:-2]  # the last two entries are special characters
                numChannels = channels.__len__()
        
        while True:
            chunk = f.read(4)
            if chunk:
                dataSample = struct.unpack('I', chunk)[0]  # data as UInt16
                data.append(dataSample)
            else:
                break
    
    # split data from file into the different channels
    # the raw data is stored such that at each time step the entires for all channels come one after the other
    # i.e. t1: d_1, d_2, ..., d_n 
    #      t2: d_1, d_2,..., d_n
    #      t3:  d_1, d_2,..., d_n
    ID = -1
    dataList = CSCsignalList()
    for i, ch in enumerate(channels):
        if not i + 1 in exclude:
            d = CSCsignal(data[i::channels.__len__()], sampleFreq)
            d.tags['channel'] = ch
            d.tags['file'] = fileName
            d.tags['path'] = os.getcwd()
        
            ID += 1
            dataList.append(ID, d)
 
    
    # add metadata to tags
    dataList.tags['file'] = fileName
    dataList.tags['path'] = os.getcwd() 
    dataList.tags['numChannels'] = numChannels
    dataList.tags['channels'] = channels
    dataList.header = header

    if showHeader:
        print header
    
    return dataList




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
    
    phidat = numpy.float_(phidat)  # make lindat to floats
    lindat = numpy.float_(lindat)  # make lindat to floats
   # lindat = lindat/lindat.max()    # normalize lindat


    # bounds of slope
    abound = [-2., 2.]

    # starting points of maximization
    Nrep = 20
    da = abs(abound[1] - abound[0]) / Nrep
    astart = min(abound)
    
    Rfunc = lambda a:-numpy.absolute(numpy.mean(numpy.exp(1j * (phidat - 2.*numpy.pi * a * lindat))))

    aopt = numpy.nan
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
    v = numpy.mean(numpy.exp(1j * (phidat - 2 * numpy.pi * aopt * lindat)))
    phi0 = numpy.angle(v)

    theta = numpy.angle(numpy.exp(2 * numpy.pi * 1j * numpy.abs(aopt) * lindat))

    thmean = numpy.angle(numpy.sum(numpy.exp(1j * theta)))
    phmean = numpy.angle(numpy.sum(numpy.exp(1j * phidat)))

    sthe = numpy.sin(theta - thmean)
    sphi = numpy.sin(phidat - phmean)

    c12 = numpy.sum(sthe * sphi)
    c11 = numpy.sum(sthe * sthe)
    c22 = numpy.sum(sphi * sphi)
    
    rho = c12 / numpy.sqrt(c11 * c22)

    lam22 = numpy.mean(sthe ** 2.*sphi ** 2)
    lam20 = numpy.mean(sphi ** 2)
    lam02 = numpy.mean(sthe ** 2)
    tval = rho * numpy.sqrt(lindat.size * lam20 * lam02 / lam22)

    p = 1 - erf(tval / numpy.sqrt(2))
    
    if display:
        fig = pl.figure(1, figsize=(10, 10))
        ax = fig.add_subplot(2, 2, 1)
        ax.plot(lindat, phidat, '.k')
        ax.plot(lindat, phidat - 2 * numpy.pi, '.k')
        ax.plot(lindat, phi0 + 2 * numpy.pi * aopt * lindat, 'r-')

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



def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len < 3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s = numpy.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]  # enlarge signal at the borders
    #    print x.shape, s.shape

    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='same')
    y = y[window_len - 1:-(window_len - 1)]  # return to original signal size
    #   print y.shape
    
    return y


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

def _read_metadata(fileName, convert2variables=False, showHeader=False):
    """
    Read the informations that may be contained in the header of
    the file.
    """
    metadata = {}
    cmd = ''
    f = open(fileName, 'r')
    for line in f.readlines():
        if line[0] == '#':
            cmd += line[1:].strip() + ';'
            if showHeader:
                print line.strip()
        elif line[0] == '%' and showHeader:  # for special comments that may be included in the header
            print line.strip()
        else:
            break
    f.close()
    exec cmd in None, metadata
    
    if showHeader:
        print ''    
    
    if convert2variables:
        for key in metadata.keys():
            exec key + '=' + str(metadata[key])

    return metadata



###################################################### CLASSES

class signal:
    
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


class CSCsignal(signal, NTsig.AnalogSignal):
    """
    A class for continuously sampled data.
    """
    
    def __init__(self, cscSignal, sampleFreq):
        
        signal.__init__(self)  # NOTE: possiblity of name clashes! AnalogSignal ahs an Attribute signal!
        self.sampleFreq = sampleFreq  # [Hz]
        dt = 1000. / sampleFreq  # [ms] expects sampling frequency in [Hz]
        
        NTsig.AnalogSignal.__init__(self, cscSignal, dt)

    def __recalc_startstop(self):
        
        self.t_start = 0.0
        self.t_stop = self.t_start + self.signal.size * self.dt
    
    def recalc_timeAxis(self):
        # just an alias to run the internal functions
        self.__recalc_startstop()
    #    self.__calcTimeAxis()
    
    def changeTimeUnit(self, newUnit='ms'):
        
        factor = _getUnitFactor(self.timeUnit, newUnit)
        
        # change the times
        self.t_start *= factor
        self.t_stop *= factor            
        self.dt *= factor
        
        self.timeUnit = newUnit
    
    
    def fft(self, display=False):
        
        n = self.signal.shape[0]
        
        if self.timeUnit == 'ms':
            binwidth = self.dt / 1000.
        elif self.timeUnit == 's':
            binwidth = self.dt
        self.sp = numpy.fft.rfft(self.signal)
        self.sp /= self.signal.size
        self.spPower = 2 * self.sp * numpy.conj(self.sp)
        self.freq = numpy.fft.fftfreq(n, d=binwidth)[0:n / 2 + 1]
        
        if display:
            self.fft_plot()
    
    
    def fft_plot(self, fig=None, ax=None):
        
        if not fig:
            fig = pl.figure(figsize=(12, 7))
        if not ax:
            ax = fig.add_subplot(111)
        
        if not hasattr(self, 'spPower'):
            self.fft()
        
        ax.plot(self.freq, self.spPower, '-')        
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power')
        ax.set_xlim(0, self.freq.max())
        
        pl.show()
        
        return fig, ax


    def filter(self, minFreq=None, maxFreq=None):
        
        if not hasattr(self, 'sp'):
            self.fft()
        
        if not minFreq:
            minFreq = self.freq.min()
        if not maxFreq:
            maxFreq = self.freq.max()
        
        # band pass filtering, zero phase
        self.sp_filtered = numpy.zeros_like(self.sp)
        for i, f in enumerate(self.freq):
            if abs(f) >= minFreq and abs(f) <= maxFreq:
                self.sp_filtered[i] = self.sp[i]
        
        # backtransform filtered signal to time space
        self.signal_filtered = numpy.fft.irfft(self.sp_filtered)
        self.signal_filtered *= self.signal.size  # rescale from normalization of the fourier components
        self.signal_filtered = self.signal_filtered[:self.signal.size]
    
    def removeFreq(self, minFreq=None, maxFreq=None):
        
        if not hasattr(self, 'sp'):
            self.fft()
        
        if not minFreq:
            minFreq = self.freq.min()
        if not maxFreq:
            maxFreq = self.freq.max()
        
        # band pass filtering, zero phase
        self.sp_purged = self.sp.copy()
        for i, f in enumerate(self.freq):
            if abs(f) >= minFreq and abs(f) <= maxFreq:
                self.sp_purged[i] = 0.0
        
        # backtransform filtered signal to time space
        self.signal_purged = numpy.fft.irfft(self.sp_purged)
        self.signal_purged *= self.signal.size  # rescale from normalization of the fourier components
        self.signal_purged = self.signal_purged[:self.signal.size]
    
    def removeMean(self):
        self.signal -= self.signal.mean()
    
    
    def spectrogram(self, minFreq=None, maxFreq=None, windowSize=4096, overlap=None, display=False):
        
        timeUnit = self.timeUnit
        if self.timeUnit != 'ms':
            self.changeTimeUnit('ms')
            
        # parameters of the spectrogram
        if not overlap:
            overlap = int(windowSize * .9)
        samplingFreq = 1. / self.dt
        
        # use matplotlibs specgram function for the spectrogram
        Pxx, freqs, t = matplotlib.mlab.specgram(self.signal, \
            NFFT=windowSize, noverlap=overlap, pad_to=windowSize * 3, Fs=samplingFreq)           
        
        freqs *= 1000  # in [Hz]
        t, freqs = numpy.meshgrid(t, freqs)
        
        if not minFreq:
            minFreq = freqs.min()
        if not maxFreq:
            maxFreq = freqs.max()
        
        indexstart = numpy.where(freqs >= minFreq)[0][0]
        indexend = numpy.where(freqs <= maxFreq)[0][-1]
        
        t = t[indexstart:indexend]
        freqs = freqs[indexstart:indexend]
        Pxx = Pxx[indexstart:indexend]
        
        if display:            
            fig = pl.figure(102, figsize=(10, 7))
            ax = fig.add_subplot(111)
            
            ax.pcolormesh(t, freqs, Pxx)
            ax.set_xlim(t.min(), t.max())
            ax.set_ylim(freqs.min(), freqs.max())
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_ylabel('Frequency [Hz]')
            
            pl.show()
        
        # switch back time unit if necessary
        if self.timeUnit != timeUnit:
            self.changeTimeUnit(timeUnit)
        
        return Pxx, freqs, t
    
    
    def time_slice(self, t_start, t_stop):
        """ Slice NeuralynxCSC between t_start and t_stop
        
        Parameters:
            t_start - begining of the new NeuralynxCSCList, in ms.
            t_stop  - end of the new NeuralynxCSCList, in ms.
        """

        assert self.t_start <= t_start
        assert self.t_stop >= t_stop
        
        t = self.time_axis()
        index1 = findNearest(t, t_start)[0]
        index2 = findNearest(t, t_stop)[0]
        self.signal = self.signal[index1:index2]
        self.recalc_timeAxis()
 
   
class CSCsignalList(signal, NTsig.AnalogSignalList):
    
    def __init__(self):
        signal.__init__(self)
        NTsig.AnalogSignalList.__init__(self, [], [], 1.)
    
    def __recalc_startstop(self):
       
        try: self.id_list
        except AttributeError:
            logging.warning("id_list is empty")
            self.signal_length = 0
        else:           
            signals = self.analog_signals.values()
            self.signal_length = len(signals[0])
            for signal in signals[1:]:
                if len(signal) != self.signal_length:
                    raise Exception("Signals must all be the same length %d != %d" % (self.signal_length, len(signal)))
        
        for csc in self:
            csc.recalc_timeAxis()
        
        self.t_start = 0.0
        self.t_stop = self.t_start + self.signal_length * csc.dt
    
    def changeTimeUnit(self, newUnit='ms'):
    
        factor = _getUnitFactor(self.timeUnit, newUnit)
        
        for csc in self:
            csc.changeTimeUnit(newUnit)
        
        self.t_start *= factor
        self.t_stop *= factor
        self.dt *= factor
        self.timeUnit = newUnit
    
    def getLimits(self):
   
        self.min = 1e100
        self.max = -1e-100
        for csc in self:
            self.min = min(self.min, csc.signal.min())
            self.max = max(self.max, csc.signal.max())
    
    def plot(self):
        
        try: self.min
        except AttributeError:
            self.getLimits()
        
        fig = pl.figure(150, figsize=(12, 7))


        ylabels = []
        last_id = self.id_list()[-1] + 1
        ax = fig.add_subplot(111)
        for id, csc in enumerate(self):
            s = csc.signal.copy()
            s -= s.mean()
            s /= self.max * .8
            ax.plot(csc.time_axis(), s + id, '-', linewidth=1)
            ylabels.append(csc.tags['channel'].split('.')[0])
        custom_plot.huebschMachen(ax)
        ax.set_xlabel('Time [' + self.timeUnit + ']')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)            
        ax.set_yticklabels(ylabels, rotation=0)
        ax.set_xlim(self.t_start, self.t_stop)
        ax.set_ylim(-1, id + 1)
        
        return fig, ax
    
    
    def fft_plot(self, freq_min=0., freq_max=60.):
        
        try: self.min
        except AttributeError:
            self.getLimits()
        
        fig = pl.figure(160, figsize=(12, 7))

        ylabels = []
        last_id = self.id_list()[-1] + 1
        ax = fig.add_subplot(111)
        for id, csc in enumerate(self):
            if not hasattr(csc, 'spPower'):
                csc.fft()
            s = csc.spPower
            s -= s.mean()
            s /= s.max() * 1.1
            ax.plot(csc.freq, s + id, '-', linewidth=1)
            ylabels.append(csc.tags['channel'].split('.')[0])
        custom_plot.huebschMachen(ax)
        ax.set_xlabel('Frequency [Hz]')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)            
        ax.set_yticklabels(ylabels, rotation=0)
        ax.set_xlim(freq_min, freq_max)
        ax.set_ylim(0, id + 1)
        
        return fig, ax

    def removeMean(self):
        for id in self.analog_signals:
            self[id].removeMean()

    def showTags(self):
        signal.showTags(self)
        for id in self.analog_signals:
            self[id].showTags()

# #    def append(self, id, signal):
# #        """
# #        Add an NeuralynxCSC object to the NeuralynxCSCList
# #
# #        Inputs:
# #            id        - the id of the new channel
# #            signal - the NeuralynxCSC object representing the new cell
# #
# #        The NeuralynxCSC object is sliced according to the t_start and t_stop times
# #        of the NeuralynxCSCList object
# #        """
# #        
# #        assert isinstance(signal, NeuralynxCSC), "An NeuralynxCSCList object can only contain NeuralynxCSC objects"
# #        if id in self.id_list():
# #            raise Exception("ID already present in NeuralynxCSCList.")
# #        else:
# #            self.signals[id] = signal

    def time_slice(self, t_start, t_stop):
        """ Slice NeuralynxCSCList between t_start and t_stop
        
        Parameters:
            t_start - begining of the new NeuralynxCSCList, in ms.
            t_stop  - end of the new NeuralynxCSCList, in ms.
        """

        for csc in self:
            csc.time_slice(t_start, t_stop)

        self.__recalc_startstop()
    
class NeuralynxCSC(CSCsignal):
    """
    A class for Neuralynx continuously sampled channel (CSC) recorded data.
    """
    
    def __init__(self, times, cscSignal, sampleFreq):
        
        signal.__init__(self)  # NOTE: possiblity of name clashes! AnalogSignal ahs an Attribute signal!
        self.times = numpy.array(times)  # time stamps provided in Neuralynx CSC file
        
        self.sampleFreq = sampleFreq  # [Hz]
        dt = 1000. / sampleFreq  # [ms] expects sampling frequency in [Hz]
        
        NTsig.AnalogSignal.__init__(self, cscSignal, dt)
        
        self.recalc_timeAxis()  

    
    def __recalc_startstop(self):
       
        self.times_start = self.times[0]
        
        try: self.tags['validSamples']
        except KeyError:
            add = 0.
        else:
            add = (self.tags['validSamples'] - 1) * self.dt
        self.times_stop = self.times[-1] + add
        
        self.t_start = 0.0
        self.t_stop = self.t_start + self.signal.size * self.dt
    
    
    def __calcTimeAxis(self):
        
        try: self.tags['validSamples']
        except KeyError:
            self.timeAxis = self.time_axis() + self.times_start
        else:
            self.timeAxis = numpy.array([])
            for i, t in enumerate(self.times[:-1]):
                self.timeAxis = numpy.append(self.timeAxis, numpy.linspace(t, self.times[i + 1], self.tags['validSamples'], endpoint=False))
            self.timeAxis = numpy.append(self.timeAxis, numpy.arange(self.times[-1], self.times[-1] + self.dt * self.tags['validSamples'], self.dt))
            # self.timeAxis = numpy.array(self.timeAxis)
            # self.timeAxis = self.timeAxis.flatten()

        # NOTES regarding time axes, t_start and t_stop etc.
        # the time stamps from the Neuralynx' CSC file might not be linear!
        #   -   therefore self.times_start and self.times_stop are place holders corresponding to
        #       the first and last times with respect to self.times, i.e., the times stamps
        #   -   t_start is put to time zero and t_stop is infered from the data, 
        #       i.e., self.t_stop = self.t_start + len(self.signal)*self.dt.
        #       They can't be offset!
        #   -   time_axis() is the time axis provided by NTsig.AnalogSignal

    
    def recalc_timeAxis(self):
        # just an alias to run the internal functions
        self.__recalc_startstop()
        self.__calcTimeAxis()
    
    
    def changeTimeUnit(self, newUnit='ms'):
        
        factor = _getUnitFactor(self.timeUnit, newUnit)
        
        # change the times
        self.t_start *= factor
        self.t_stop *= factor            
        self.dt *= factor
        
        self.times_start *= factor
        self.times_stop *= factor
        self.times *= factor
        self.timeAxis *= factor       
        
        self.timeUnit = newUnit

    
    def hilbertTransform(self, display=False, justPhase=True):
        """ Compute the hilbert transform and from that the analytical signal 
        of the filtered analog signal and return the phase and 
        absolute values of the later.
        
        Parameters:
            display ... optional display of a figure with the spiketrains and ints convolved counterpart
            justPhase ... just return the hilbert phase"""
    
        try: self.signal_filtered
        except AttributeError:
            self.filter()
            print "NOTE: Analog signal was not filtered before. Filtered it now using default values."
        
        self.analyticalTrain = scipy.signal.hilbert(self.signal_filtered.real)
        self.hilbertAbsolute = numpy.absolute(self.analyticalTrain)
        self.hilbertPhase = numpy.angle(self.analyticalTrain, deg=True)
        for i in range(self.hilbertPhase.shape[-1]):
            if self.hilbertPhase[i] < 0:
                self.hilbertPhase[i] += 360
        
        if display:
            fig = pl.figure(103, figsize=(12, 7))
            ax = fig.add_subplot(311)   

            ax.plot(self.timeAxis, self.analyticalTrain.real)
            ax.plot(self.timeAxis, self.analyticalTrain.imag)
            ax.set_xlim(self.t_start, self.t_stop)
            ax.legend(['Signal', 'Hilbert transform'])

            ax = fig.add_subplot(312)
            ax2 = ax.twinx()
            ax.plot(self.timeAxis, self.hilbertAbsolute)
            ax2.plot(self.timeAxis, self.hilbertPhase, 'g')
            ax.set_xlim(self.t_start, self.t_stop) 
            ax2.set_ylim(0, 360)
            ax2.set_yticks(range(0, 410, 90))
            ax2.set_yticklabels(range(0, 410, 90))
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_ylabel('Absolute value')
            ax2.set_ylabel('Phase [deg]')

            ax = fig.add_subplot(313)
            ax2 = ax.twinx()
            ax.plot(self.timeAxis, self.analyticalTrain.real)
            ax2.plot(self.timeAxis, self.hilbertPhase, 'g')
            ax.set_xlim(self.t_start, self.t_stop) 
            ax2.set_ylim(0, 360) 
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_ylabel('Signal')
            ax2.set_ylabel('Phase [deg]')

            pl.show()
        
        if justPhase:
            return self.hilbertPhase
        else:
            return self.hilbertPhase, self.hilbertAbsolute, self.analyticalTrain  

    
    def plot(self, fig=None, ax=None):
        
        if not fig:
            fig = pl.figure(figsize=(12, 7))
        if not ax:
            ax = fig.add_subplot(111)
        
        ax.plot(self.time_axis(), self.signal, '-')
        
        ax.set_xlabel('Time [' + self.timeUnit + ']')
        ax.set_xlim(self.t_start, self.t_stop)
        
        pl.show()
        
        return fig, ax
    
    
    def purge(self):
        """ Remove values that fill up the possible range.
        """
        
        ds = numpy.diff(self.signal)
        indices, = numpy.where(ds == 0)
        indices = numpy.concatenate((indices - 1, indices, indices + 1))
        self.signal = numpy.delete(self.signal, indices)
        
        # to do properly:
        # self.times ... purge as well

        self.recalc_timeAxis()
    
    def resample(self, new_dt):
        if new_dt > self.dt:
            bins = int(new_dt / self.dt)  # NOTE: Not interpolating!!!!!!!
            print bins
            print new_dt / self.dt
            self.signal = self.signal[::bins]
            self.dt = new_dt
        else:
            print "Not implemented yet."
    
    
    def time_offset(self, offset):
        """
        Add a time offset to the NeuralynxCSC object.
        """
        self.times += offset
        self.times_start += offset
        self.times_stop += offset

    def time_slice(self, t_start, t_stop):
        """ Slice NeuralynxCSC between t_start and t_stop
        
        Parameters:
            t_start - begining of the new NeuralynxCSCList, in ms.
            t_stop  - end of the new NeuralynxCSCList, in ms.
        """

        index1 = numpy.searchsorted(self.times, t_start)
        index2 = numpy.searchsorted(self.times, t_stop)
        if index2 > self.times.shape[0]:
            index2 = self.times.shape[0]
        self.times = self.times[index1:index2]
        self.signal = self.signal[index1 * self.tags['validSamples']:index2 * self.tags['validSamples']]
        self.recalc_timeAxis()
    

    
    
class NeuralynxCSCList(CSCsignalList):
    
    def __init__(self):
        CSCsignalList.__init__(self)
    
    def __recalc_startstop(self):
       
        try: self.id_list
        except AttributeError:
            logging.warning("id_list is empty")
            self.signal_length = 0
        else:
            start_times = numpy.array([self.analog_signals[id].timeAxis[0] for id in self.id_list()], numpy.float32)
            stop_times = numpy.array([self.analog_signals[idx].timeAxis[-1] for idx in self.id_list()], numpy.float32)
            
            signals = self.analog_signals.values()
            self.signal_length = len(signals[0])
            for signal in signals[1:]:
                if len(signal) != self.signal_length:
                    raise Exception("Signals must all be the same length %d != %d" % (self.signal_length, len(signal)))
        
        for csc in self:
            csc.__recalc_startstop()
        
        self.t_start = 0.0
        self.t_stop = self.t_start + self.signal_length * csc.dt
    
    
    def plot(self):
        
        try: self.min
        except AttributeError:
            self.getLimits()
        
        fig = pl.figure(150, figsize=(12, 7))


        ylabels = []
        last_id = self.id_list()[-1] + 1
        ax = fig.add_subplot(111)
        for id, csc in enumerate(self):
            s = csc.signal.copy()
            s -= s.mean()
            s /= self.max * .8
            ax.plot(csc.time_axis(), s + id, '-', linewidth=1)
            ylabels.append(csc.tags['file'].split('.')[0])
        custom_plot.huebschMachen(ax)
        ax.set_xlabel('Time [' + self.timeUnit + ']')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)            
        ax.set_yticklabels(ylabels, rotation=0)
        ax.set_xlim(self.t_start, self.t_stop)
        ax.set_ylim(-1, id + 1)
        
        return fig, ax
    
    
    def fft_plot(self, freq_min=0., freq_max=60.):
        
        try: self.min
        except AttributeError:
            self.getLimits()
        
        fig = pl.figure(160, figsize=(12, 7))

        ylabels = []
        last_id = self.id_list()[-1] + 1
        ax = fig.add_subplot(111)
        for id, csc in enumerate(self):
            if not hasattr(csc, 'spPower'):
                csc.fft()
            s = csc.spPower
            s -= s.mean()
            s /= s.max() * 1.1
            ax.plot(csc.freq, s + id, '-', linewidth=1)
            ylabels.append(csc.tags['file'].split('.')[0])
        custom_plot.huebschMachen(ax)
        ax.set_xlabel('Frequency [Hz]')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)            
        ax.set_yticklabels(ylabels, rotation=0)
        ax.set_xlim(freq_min, freq_max)
        ax.set_ylim(0, id + 1)
        
        return fig, ax



# #    def append(self, id, signal):
# #        """
# #        Add an NeuralynxCSC object to the NeuralynxCSCList
# #
# #        Inputs:
# #            id        - the id of the new channel
# #            signal - the NeuralynxCSC object representing the new cell
# #
# #        The NeuralynxCSC object is sliced according to the t_start and t_stop times
# #        of the NeuralynxCSCList object
# #        """
# #        
# #        assert isinstance(signal, NeuralynxCSC), "An NeuralynxCSCList object can only contain NeuralynxCSC objects"
# #        if id in self.id_list():
# #            raise Exception("ID already present in NeuralynxCSCList.")
# #        else:
# #            self.signals[id] = signal



class NeuralynxEvents(signal):
    """
    A class for Neuralynx event data.
    """
    
    def __init__(self, times, eventStrings):
        signal.__init__(self)
        self.times = numpy.array(times)
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

 

class inputSignals(signal):
    
    types = {0:'exi', 1:'inh', 2:'el_soma', 3:'condel_soma', 4:'noise_soma', 5:'el_dend', 6:'condel_dend', 7:'noise_dend'}
    
    def __init__(self, fileName):
        signal.__init__(self)
        self.inputParams = _read_metadata(fileName)
        for key in self.inputParams.keys():  # convert dictionary to variables of the inputSignals object
            exec 'self.' + key + '=' + str(self.inputParams[key])
            
        self.fileName = fileName
        self.inputs = []
        self.loadInputs(fileName)
        self.numInputs = self.inputs.__len__()
        self.timeOffset = 0.0


    def loadInputs(self, fileName):  
        data = numpy.loadtxt(fileName)
        
        for i in range(self.inputParams['first_id'], self.inputParams['last_id'] + 1):
            indices = numpy.nonzero(data[:, 2] == i)
            self.inputs.append(data[indices[0], 0:2])
    
    def plotInputs(self):
        fig = pl.figure(figsize=(10, 10))
        for i in range(self.numInputs):
            ax = fig.add_subplot(self.numInputs + 1, 1, i + 1)
            if self.inputs[i].size:
                ax.plot(self.inputs[i][:, 0], self.inputs[i][:, 1])
            
            ax.set_xlim(self.getMinTime(), self.getMaxTime())
            ax.set_ylabel(self.types[i])
            if i < self.numInputs - 1:
                ax.set_xticklabels([])
        ax = fig.add_subplot(self.numInputs + 1, 1, self.numInputs + 1)
        ax.plot(self.summedInputs()[:, 0], self.summedInputs()[:, 1])
        ax.set_xlim(self.getMinTime(), self.getMaxTime())
        ax.set_ylabel('sum')
        
        ax.set_xlabel('Time [' + self.timeUnit + ']')
        pl.show()
    
    def offsetTime(self, offset):
        self.timeOffset = offset
        self.inputs[:, 0] += offset
    
    def getMinTime(self):
        minTime = numpy.inf
        for i in self.inputs:
            if i.size and minTime > i[:, 0].min():
                minTime = i[:, 0].min()
        return minTime
    
    def getMaxTime(self):
        maxTime = -numpy.inf
        for i in self.inputs:
            if i.size and maxTime < i[:, 0].max():
                maxTime = i[:, 0].max()
        return maxTime
    
    def summedInputs(self):
        summe = self.inputs[0].copy()  
        for index0 in range(1, self.inputs.__len__()):
            inp = self.inputs[index0]
            for j, time in enumerate(inp[:, 0]):
                index = summe[:, 0].searchsorted(time)
                if summe.shape[0] > index and summe[index, 0] == time:
                    summe[index, 1] += inp[j, 1]   
                else:
                    summe = numpy.insert(summe, index, inp[j], axis=0)
        return summe

class placeCell_spikezug(signal, NTsig.SpikeTrain):
    """ Spike train class for a place modulated cell.
    
    It extends the NeuroTools SpikeTrain class.
    """
    
    def __init__(self, spike_times, t_start=None, t_stop=None, timeUnit='ms'):
        signal.__init__(self)
        NTsig.SpikeTrain.__init__(self, spike_times, t_start=None, t_stop=None)

    def getLeftAndRightwardSpikes(self, traj, onlyRunning=False):
        
        inclPhases = True
        if not hasattr(self , 'spike_phases'):
            print 'WARNING: Calculating without spike phases! If spike phases are needed \
                caluculate spike phases have to be calculated first!'
            inclPhases = False
        
        if onlyRunning:
            rechts_spikeTimes = numpy.copy(self.run_spikeTimes)
            rechts_spikePlaces = numpy.copy(self.run_spikePlaces)
            links_spikeTimes = numpy.copy(self.run_spikeTimes)
            links_spikePlaces = numpy.copy(self.run_spikePlaces)
            
            if inclPhases:
                rechts_spikePhases = numpy.copy(self.run_spikePhases)
                links_spikePhases = numpy.copy(self.run_spikePhases)
        else:
            rechts_spikeTimes = numpy.copy(self.spike_times)
            rechts_spikePlaces = numpy.copy(self.spikePlaces)
            links_spikeTimes = numpy.copy(self.spike_times)
            links_spikePlaces = numpy.copy(self.spikePlaces)
            
            if inclPhases:
                rechts_spikePhases = numpy.copy(self.spike_phases)
                links_spikePhases = numpy.copy(self.spike_phases)
        
        for j, time in enumerate(self.run_spikeTimes):
            if numpy.ma.count_masked(numpy.ma.masked_values(traj.rechts_times, time, atol=traj.dt)) == 0:
                rechts_spikeTimes[j] = numpy.nan
                rechts_spikePlaces[j] = numpy.nan
                if inclPhases:
                    rechts_spikePhases[j] = numpy.nan
            if numpy.ma.count_masked(numpy.ma.masked_values(traj.links_times, time, atol=traj.dt)) == 0:
                links_spikeTimes[j] = numpy.nan
                links_spikePlaces[j] = numpy.nan
                if inclPhases:
                    links_spikePhases[j] = numpy.nan


        # remove nan's from Spike-arrays!!
        if rechts_spikeTimes.size > 0:  # sonst fehler, wenn array empty
            self.rechts_spikeTimes = rechts_spikeTimes[~numpy.isnan(rechts_spikeTimes)]
            self.rechts_spikePlaces = rechts_spikePlaces[~numpy.isnan(rechts_spikePlaces).any(1)]
            if inclPhases:
                self.rechts_spikePhases = rechts_spikePhases[~numpy.isnan(rechts_spikePhases)]
        if links_spikeTimes.size > 0:  # sonst fehler, wenn array empty
            self.links_spikeTimes = links_spikeTimes[~numpy.isnan(links_spikeTimes)]
            self.links_spikePlaces = links_spikePlaces[~numpy.isnan(links_spikePlaces).any(1)]
            if inclPhases:
                self.links_spikePhases = links_spikePhases[~numpy.isnan(links_spikePhases)]

    
    def getSpikePlaces(self, traj, interp=True):
        """ Return the place at which a certain spike occured.
        """
        
        if interp:
            print "NOTE: Using linear interpolation to get the places."
        
        spikePlaces = []
        for spike in self.spike_times:
            if spike > traj.times[-1]:
                print 'WARNING: at time', spike, 'there was a spike after the trajectory was over!'
                spikePlaces.append(numpy.ones(3) * (0))
            else:
                spikePlaces.extend(traj.getPlaceFromTime(spike, interp=interp))
        
        self.spikePlaces = numpy.array(spikePlaces)
    
    
    def getRunningSpikes(self, traj, threshspeed=0.01):
        """ Get spikes that occured above certain running speed.
        """
        
        inclPhases = True
        if not hasattr(self , 'spike_phases'):
            print 'WARNING: Calculating without spike phases! If spike phases are needed caluculate them before!'
            inclPhases = False
        
        inclHeadDir = True
        if not hasattr(self , 'spike_headDirections'):
            print 'WARNING: Calculating without spike head directions! If spike head directions are needed caluculate them before!'
            inclHeadDir = False
        
        if inclPhases:
            self.run_spikePhases = numpy.copy(self.spike_phases)
        if inclHeadDir:
            self.run_spikeHeadDirections = numpy.copy(self.spike_headDirections)
        
        
        self.run_spikeTimes = numpy.copy(self.spike_times)
        if not hasattr(self , 'spikePlaces'):
            self.getSpikePlaces(traj)
        self.run_spikePlaces = numpy.copy(self.spikePlaces)
        
        if not hasattr(traj , 'run_times'):
            print 'NOTE: Now calculating running traj >', threshspeed, traj.spaceUnit + '/' + traj.timeUnit 
            traj.getRunningTraj(threshspeed=threshspeed)
        
        for j, time in enumerate(self.spike_times):
            if numpy.ma.count_masked(numpy.ma.masked_values(traj.run_times, time, atol=traj.dt)) == 0:
                self.run_spikeTimes[j] = numpy.nan
                self.run_spikePlaces[j] = numpy.nan
                if inclPhases:
                    self.run_spikePhases[j] = numpy.nan
                if inclHeadDir:
                    self.run_spikeHeadDirections[j] = numpy.nan

        # remove nans from arrays!
        if self.run_spikeTimes.size > 0:  # sonst fehler, wenn array empty
            self.run_spikeTimes = self.run_spikeTimes[~numpy.isnan(self.run_spikeTimes)]
            self.run_spikePlaces = self.run_spikePlaces[~numpy.isnan(self.run_spikePlaces).any(1)]
            if inclPhases:
                self.run_spikePhases = self.run_spikePhases[~numpy.isnan(self.run_spikePhases)]
            if inclHeadDir:
                self.run_spikeHeadDirections = self.run_spikeHeadDirections[~numpy.isnan(self.run_spikeHeadDirections)]
    
    
    def getSpikePhases(self, lfp):
        """ Get firing phase with respect to LFP.
        """
        
        spikes = self.spike_times
        
        if not hasattr(lfp , 'hilbertPhase'):
            lfp.hilbertTransform()
        
        phases = []
        for spike in spikes:
            # phases.append(lfp.hilbertPhase[(lfp.timeAxis).searchsorted(spike)])
            phases.append(lfp.hilbertPhase[findNearest(lfp.timeAxis, spike)[0]])
        self.spike_phases = numpy.array(phases)
        
        return self.spike_phases
        
        
    def getSpikeHeadDirection(self, traj, interp=True):
        """ Return the head direction at which a certain spike occured.
        """
        
        if interp:
            print "NOTE: Using interpolation to get head directions."
        
        spike_headDirections = []
        for spike in self.spike_times:
            if spike > traj.times[-1]:
                print 'WARNING: at time', spike, 'there was a spike after the trajectory was over!'
                spike_headDirections.append(0)
            else:
                spike_headDirections.extend(traj.getHeadDirectionFromTime(spike, interp=interp))
    
        self.spike_headDirections = numpy.array(spike_headDirections)


    def plotSpikesvsTraj(self, traj, onlyRunning=False, showHeadDir=False, fig=None, ax=None):
        """ Get firing field for pseudo color plotting.
        """
        
        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)
        
        if not hasattr(self , 'spikePlaces'):
            self.getSpikePlaces(traj)
    
        ax.plot(traj.places[:, 0], traj.places[:, 1], linewidth=1, color=numpy.ones(3) * .25)  # original/full trajectory 
        
        if showHeadDir:
            hd = traj.getHeadDirection()
            ax.quiver(traj.places[:, 0], traj.places[:, 1], numpy.cos(hd), numpy.sin(hd), pivot='mid', color='g', units='y')


        if not onlyRunning:
            ax.plot(self.spikePlaces[:, 0], self.spikePlaces[:, 1], color='#FFA1A1', marker='.', linestyle='None')  # , label='all spikes')  # all spiketimes
        ax.plot(self.run_spikePlaces[:, 0], self.run_spikePlaces[:, 1], 'r.')  # , label='only running spikes') # only running-spikes

        if showHeadDir:
            if not onlyRunning:
                ax.quiver(self.spikePlaces[:, 0], self.spikePlaces[:, 1], \
                    numpy.cos(self.spike_headDirections), \
                    numpy.sin(self.spike_headDirections), pivot='mid', color=[.1, .1, 1], units='dots')
            ax.quiver(self.run_spikePlaces[:, 0], self.run_spikePlaces[:, 1], \
                numpy.cos(self.run_spikeHeadDirections), \
                numpy.sin(self.run_spikeHeadDirections), pivot='mid', color=[0, 0, 1], units='dots')

        ax.set_xlim(traj.xlim + numpy.diff(traj.xlim) * .05 * [-1, 1])
        ax.set_ylim(traj.ylim + numpy.diff(traj.ylim) * .1 * [-1, 1])
    
        ax.set_xlabel('x position (' + traj.spaceUnit + ')')
        ax.set_ylabel('y position (' + traj.spaceUnit + ')')
        
        return fig, ax
    
    
    def plotField(self, traj, bin_size=.05, onlyRunning=False, fig=None, ax=None):
        """ Get firing field for pseudo color plotting.
        """
        
        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)
        
        if onlyRunning:
            spikes = self.run_spikePlaces[:, :2]
        else:
            spikes = self.spikePlaces[:, :2]
        
        # bin track!
        track_x1 = numpy.arange(traj.xlim[0], traj.xlim[1], bin_size)
        track_x2 = numpy.arange(traj.xlim[0] + bin_size, traj.xlim[1] + bin_size, bin_size)
        track_y1 = numpy.arange(traj.ylim[0], traj.ylim[1], bin_size)
        track_y2 = numpy.arange(traj.ylim[0] + bin_size, traj.ylim[1] + bin_size, bin_size)

        # create polygones & count spikes
        spike_number = []
        for l, ySP in enumerate(track_y1):  #
            for j, xSP in enumerate(track_x1):
                pol = [[track_x1[j], track_y1[l]], [track_x2[j], track_y1[l]], \
                    [track_x1[j], track_y2[l]], [track_x2[j], track_y2[l]]]  # polygon erzeugen
                # pol.append([[track_x1[j],track_y1[l]], [track_x2[j],track_y1[l]], [track_x1[j],track_y2[l]], [track_x2[j],track_y2[l]]])
                spike_number.append(numpy.sum(mpl_nxutils.points_inside_poly(spikes, pol)))  # count number of spikes in polygon
        X, Y = numpy.meshgrid(track_x1, track_y1)
        spike_number = numpy.array(spike_number).reshape(len(track_y1), len(track_x1))  # reshape & array spike_number list
        # X,Y = numpy.meshgrid(numpy.arange(traj.xlim[0],traj.xlim[1]+bin_size, bin_size), numpy.arange(traj.ylim[0], traj.ylim[1]+bin_size, bin_size))


        ax.pcolor(X, Y, spike_number)

        ax.set_xlim(traj.xlim[0], traj.xlim[1] - bin_size)
        ax.set_ylim(traj.ylim[0], traj.ylim[1] - bin_size)

        ax.set_xlabel('x position (' + traj.spaceUnit + ')')
        ax.set_ylabel('y position (' + traj.spaceUnit + ')')
        
        return fig, ax
    
    
    def time_slice(self, t_start, t_stop):

        if hasattr(t_start, '__len__'):
            if len(t_start) != len(t_stop):
                raise ValueError("t_start has %d values and t_stop %d. They must be of the same length." % (len(t_start), len(t_stop)))
            mask = False
            for t0, t1 in zip(t_start, t_stop):
                mask = mask | ((self.spike_times >= t0) & (self.spike_times <= t1))
            t_start = t_start[0]
            t_stop = t_stop[-1]
        else:
            mask = (self.spike_times >= t_start) & (self.spike_times <= t_stop)
        spikes = numpy.extract(mask, self.spike_times)
        
        return placeCell_spikezug(spikes, t_start, t_stop)



class spikezug(signal, NTsig.SpikeTrain):
    """
    Spike train class that extends the NeuroTools SpikeTrain class.
    """
    
    def __init__(self, spike_times, timeUnit='ms'):
        signal.__init__(self)
        NTsig.SpikeTrain.__init__(self, spike_times)

    def getSpikePhasesFreq(self, osciFreq=1, osciOffset=0):
        """ calculate spike phases with an fixed oscillation frequency """
        if self.timeUnit == 'ms':
            T = 1 / osciFreq * 1000  # period length [ms]
        elif self.timeUnit == 's':
            T = 1 / osciFreq  # period length [s]
        else:
            print 'WARNING: Time unit not adequately supported.'
        
        self.spike_phases = numpy.zeros(self.spike_times.shape[0])
        for i, spike in enumerate(self.spike_times):
            self.spike_phases[i] = spike / T
            self.spike_phases[i] -= int(spike / T)
            self.spike_phases[i] *= 360
        return self.spike_phases
    
    def getSpikePhasesT(self, t_start=0, t_stop=0):
        """ calculate spike phases for an cycle lasting from t_start to t_stop """
        T = t_stop - t_start
        timeSlicedTrain = self.time_slice(t_start, t_stop).spike_times
        spike_phases = numpy.zeros(timeSlicedTrain.shape[0])
        for i, spike in enumerate(timeSlicedTrain):
            spike_phases[i] = (spike - t_start) / T
            spike_phases[i] *= 360
        return spike_phases


class spikezugList(signal, NTsig.SpikeList):
    """
    A class that extends the NeuroTools SpikeList class.
    """
    
    def __init__(self, fileName=None, timeUnit='ms', **kwargs):
        signal.__init__(self)
        self.timeUnit = timeUnit
        if fileName:
            data = numpy.loadtxt(fileName)
            data = numpy.hstack((data[:, [1]], data[:, [0]]))  # swap left and right columns in data, necessary to use SpikeList.__init__ properly
            data = [tuple(t) for t in data]
            self.meta = _read_metadata(fileName)
            NTsig.SpikeList.__init__(self, data, range(self.meta['first_id'], self.meta['last_id'] + 1))
        else:
            NTsig.SpikeList.__init__(self, [], [], t_start=kwargs['t_start'], \
                t_stop=kwargs['t_stop'], dims=kwargs['dims'])
    
    def __calcTimeAxis(self):
        """
        Calculate the time axis for the LFP, hilbert phases etc.
        """
        self.timeAxis = self.time_axis(self.binwidth)[:-1]
        # self.timeAxis = self.time_axis(binwidth)
        
    def __recalc_startstop(self):
        start_times = numpy.array([self.spiketrains[id].spike_times[0] for id in self.id_list], numpy.float32)
        self.t_start = numpy.min(start_times)
        for id in self.spiketrains.keys():
            self.spiketrains[id].t_start = self.t_start
        
        stop_times = numpy.array([self.spiketrains[idx].spike_times[-1] for idx in self.id_list], numpy.float32)
        self.t_stop = numpy.max(stop_times)
        for id in self.spiketrains.keys():
            self.spiketrains[id].t_stop = self.t_stop


    def addTags(self, id, **kwargs):
        """
        For adding and changing comments/tags to the spiketrains.
        """
        self.tags[id] = {}
        for item in kwargs:
            self.tags[id][item] = kwargs[item]

    def showTags(self):
        """
        For showing the comments/tags to the spiketrains on the command line.
        """
        for id in self.tags:
            print 'cell #', id
            for item in self.tags[id]:
                print ' ', item, ':', self.tags[id][item]
            print ''
 
    
    def changeTimeUnit(self, newUnit='ms'):
        
        factor = _getUnitFactor(self.timeUnit, newUnit)
        
        if factor != 1.0:
        
            # change the times
            for id in self.spiketrains:
                self.spiketrains[id].spike_times *= factor
                self.spiketrains[id].t_start *= factor
                self.spiketrains[id].t_stop *= factor
                self.spiketrains[id].timeUnit = newUnit          
            self.__recalc_startstop()
            
            try: self.timeAxis
            except AttributeError:
                pass
            else:
                self.timeAxis *= factor
            
            try: self.binwidth                             
            except AttributeError:
                pass
            else:
                self.binwidth *= factor
                
            self.timeUnit = newUnit

    
    def convolveSpikes(self, kernel, kernel_offset=0, binwidth=1., display=False):
        """Convolve the spiketrains in the spikezugList with a given kernel
        
        Parameters:
            kernel_offset ... offset the kernel by a certain number of bins
            binwidth ... time between bins, NOTE: also affects the convolution, since the kernel is just defined in bins!
            display ... optional display of a figure with the spiketrains and ints convolved counterpart
        """
        self.convolvedSpikeTrains = []
        spikeTrainHist = []
        if inspect.ismethod(self.id_list):
            IDlist = self.id_list()
        else:
            IDlist = self.id_list
        
# #        try: self.timeAxis
# #        except AttributeError:
# #            self.__calcTimeAxis()
# #        
# #        try: self.binwidth                              # binwidth
# #        except AttributeError:
# #            self.binwidth = binwidth
# #        else:
# #            if self.binwidth != binwidth:
# #                self.binwidth = binwidth
# #                self.__calcTimeAxis()

        self.binwidth = binwidth
        self.__calcTimeAxis()  # time axis for the LFP, hilbert phases etc.
            
        for id in IDlist:
            train = self.__getitem__(id)
# #  NOTE: lines below moved upwards and using the SpikeList's not SpikeTrain's time_axis()
# #            try: self.timeAxis
# #            except AttributeError:
# #                self.timeAxis=train.time_axis(binwidth)[:-1]
            spikeTrainHist.append(train.time_histogram(self.binwidth, normalized=False))
            conv = numpy.convolve(spikeTrainHist[-1], kernel, mode='full')
            # conv /= self.timeAxis.shape[-1]            # normalize convolution value
            conv = conv[kernel_offset:-kernel.shape[-1] + kernel_offset + 1]  # reduce size so that conv has right shape and spikes and convolution are appropriately aligned
                                            # NOTE: Problem in case kernel_offset>kernel.shape[-1], needs to be rewritten!
            self.convolvedSpikeTrains.append(conv)        
        
        if display:
            fig = pl.figure(100, figsize=(12, 6))
            
            ax = fig.add_subplot(211)
            pos = ax.get_position()
            pos.x1 = .65
            ax.set_position(pos)
            for hist in spikeTrainHist:
                ax.plot(self.timeAxis, hist, color=[.2, .2, .2])
            ax.set_ylabel('Spikes')
            ax.set_xlim(self.t_start, self.t_stop)
            
            ax = fig.add_subplot(212)
            pos = ax.get_position()
            pos.x1 = .65
            ax.set_position(pos)
            for convTrain in self.convolvedSpikeTrains:
                ax.plot(self.timeAxis, convTrain)
            # ax.set_ylim(0, maxiPow*1.1)
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_xlim(self.t_start, self.t_stop)
            ax.set_ylabel('Convolved spikes')
            # ax.legend(['Raw signal','Filtered signal '+str(minFreq)+'-'+str(maxFreq)+' Hz'])
            
            ax = fig.add_subplot(111)
            ax.set_position([pos.x1 + .1, .1, .2, pos.height * 2.2])
            pl.plot((numpy.arange(0, kernel.shape[0]) - kernel_offset) * self.binwidth, kernel)
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_title('Kernel')
            
            pl.show()
    
    def sumConvolvedSpikeTrains(self, first_id=None, last_id=None, display=False):
        """Sum the convolved spiketrains.
                
        Parameters:
            first_id
            last_id
            display ... optional display of a figure with the spiketrains and ints convolved counterpart
        """
        if inspect.ismethod(self.id_list):
            IDlist = self.id_list()
        else:
            IDlist = self.id_list
        if not first_id:
            first_id = IDlist[0]
        if not last_id:
            last_id = IDlist[-1]
        self.summedConvolvedSpikeTrains = sum(self.convolvedSpikeTrains[first_id:last_id + 1])
        
        if display:
            fig = pl.figure(101, figsize=(12, 5))
            
            ax = fig.add_subplot(111)
            for convTrain in self.convolvedSpikeTrains:
                ax.plot(self.timeAxis, convTrain)
            ax.plot(self.timeAxis, self.summedConvolvedSpikeTrains, linewidth=3)
            # ax.set_ylim(0, maxiPow*1.1)
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_xlim(self.t_start, self.t_stop)
            ax.set_ylabel('Summed convolved spikes')
            # ax.legend(['Raw signal','Filtered signal '+str(minFreq)+'-'+str(maxFreq)+' Hz'])
            
            pl.show()
        
        return self.summedConvolvedSpikeTrains

    
    def filterSummedConvolvedSpikeTrains(self, minFreq=None, maxFreq=None, timeRange=None, display=False, incRes=1):
        """For filtering the sum of the convolved spiketrains.
                
        Parameters:
            minFreq, maxFreq ... band pass filter on the interval [minFreq, maxFreq]
            timeRange ... tuple (start time, end time)
            display ... optional display of a figure with the spiketrains and ints convolved counterpart
            incRes ... factor used for increasing the number of entries in self.summedConvolvedSpikeTrains,
                        i.e., increasing the resolution on the frequency axis
                        NOTE: has to be an integer!
        """
        try: self.convolvedSpikeTrains
        except AttributeError:
            self.convolveSpikes(numpy.array([1]))
            
        try: self.summedConvolvedSpikeTrains
        except AttributeError:
            self.sumConvolvedSpikeTrains()
            
        if not timeRange:
            timeRange = (self.timeAxis[0], self.timeAxis[-1])
        # convert timeRange from time to indices
        indexRange = (numpy.searchsorted(self.timeAxis, timeRange[0]), \
                numpy.searchsorted(self.timeAxis, timeRange[1]))
        self.summedConvolvedSpikeTrains = self.summedConvolvedSpikeTrains[indexRange[0]:indexRange[1] + 1]
        self.timeAxis = self.timeAxis[indexRange[0]:indexRange[1] + 1]
        
        # perform fft
        n = self.summedConvolvedSpikeTrains.shape[-1] * incRes  # trick to get better frequency resolution
        n = int(2 ** numpy.ceil(numpy.log(n) / numpy.log(2)))  # to get n = 2**x
        # print n
        self.sp = numpy.fft.rfft(self.summedConvolvedSpikeTrains, n=n)
        self.sp /= self.timeAxis.shape[-1]  # normalize fourier components to the number of sample points in signal
        self.freq = numpy.fft.fftfreq(n, d=self.binwidth / 1000.)[0:n / 2 + 1]
        
        if not minFreq:
            minFreq = self.freq.min()
        if not maxFreq:
            maxFreq = self.freq.max()
        
        # band pass filtering, zero phase
        self.sp_filtered = numpy.zeros_like(self.sp)
        for i, f in enumerate(self.freq):
            if abs(f) >= minFreq and abs(f) <= maxFreq:
                self.sp_filtered[i] = self.sp[i]
        
        # backtransform filtered signal to time space
        self.filteredSummedConvolvedSpikeTrains = numpy.fft.irfft(self.sp_filtered)
        self.filteredSummedConvolvedSpikeTrains *= self.timeAxis.shape[-1]  # rescale from normalization of the fourier components
        self.filteredSummedConvolvedSpikeTrains = self.filteredSummedConvolvedSpikeTrains.real
        self.filteredSummedConvolvedSpikeTrains = self.filteredSummedConvolvedSpikeTrains[:self.summedConvolvedSpikeTrains.shape[-1]]

        # determine power
        self.spPower = 2 * self.sp * self.sp.conj()
        self.sp_filteredPower = 2 * self.sp_filtered * self.sp_filtered.conj()
        
        # plot
        if display:
            fig = pl.figure(102, figsize=(12, 7))
            
            ax = fig.add_subplot(211)
            ax.plot(self.freq, self.spPower, '.-')  # times 2, since we are just on the real axis
            ax.plot(self.freq, self.sp_filteredPower)  # times 2, since we are just on the real axis
            ax.set_xlim(0, min(self.freq.max(), maxFreq * 10))
            maxiPow = self.spPower[1:].max()  # discard first power entry
            ax.set_ylim(0, maxiPow * 1.1)
            ax.set_xlabel('Frequency [Hz]')
            ax.set_ylabel('Power')
            ax.legend(['Raw signal', 'Filtered signal ' + str(minFreq) + '-' + str(maxFreq) + ' Hz'])
            
            ax = fig.add_subplot(212)
            ax.plot(self.timeAxis, self.summedConvolvedSpikeTrains)
            ax.plot(self.timeAxis, self.filteredSummedConvolvedSpikeTrains)
            # ax.plot(self.timeAxis, self.filteredSummedConvolvedSpikeTrains.real)
            # ax.plot(self.timeAxis, self.filteredSummedConvolvedSpikeTrains.imag)
            ax.set_xlim(self.t_start, self.t_stop) 
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.legend(['Raw signal', 'Filtered signal ' + str(minFreq) + '-' + str(maxFreq) + ' Hz'])
            pl.show()
        
        return self.filteredSummedConvolvedSpikeTrains

    
   

    def hilbertTransform(self, display=False, justPhase=True):
        """ Compute the hilbert transform and from that the analytical signal 
            of the filtered convolved spike trains and return the phase and 
            absolute values of the later.
            
            Parameters:
                display ... optional display of a figure with the spiketrains and ints convolved counterpart
                justPhase ... just return the hilbert phase"""
         
        try: self.filteredSummedConvolvedSpikeTrains
        except AttributeError:
            self.filterSummedConvolvedSpikeTrains()
        
            
        self.analyticalTrain = scipy.signal.hilbert(self.filteredSummedConvolvedSpikeTrains.real)
        self.hilbertPhase = numpy.angle(self.analyticalTrain, deg=True)
        for i in range(self.hilbertPhase.shape[-1]):
            if self.hilbertPhase[i] < 0:
                self.hilbertPhase[i] += 360
        self.hilbertAbsolute = numpy.absolute(self.analyticalTrain)
            
        if display:
            fig = pl.figure(103, figsize=(12, 7))
            ax = fig.add_subplot(311)   
            
            ax.plot(self.timeAxis, self.analyticalTrain.real)
            ax.plot(self.timeAxis, self.analyticalTrain.imag)
            ax.set_xlim(self.t_start, self.t_stop)
            ax.legend(['Signal', 'Hilbert transform'])
            
            ax = fig.add_subplot(312)
            ax2 = ax.twinx()
            ax.plot(self.timeAxis, self.hilbertAbsolute)
            ax2.plot(self.timeAxis, self.hilbertPhase, 'g')
            ax.set_xlim(self.t_start, self.t_stop) 
            ax2.set_ylim(0, 360)
            ax2.set_yticks(range(0, 410, 90))
            ax2.set_yticklabels(range(0, 410, 90))
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_ylabel('Absolute value')
            ax2.set_ylabel('Phase [deg]')
            
            ax = fig.add_subplot(313)
            ax2 = ax.twinx()
            ax.plot(self.timeAxis, self.analyticalTrain.real)
            ax2.plot(self.timeAxis, self.hilbertPhase, 'g')
            ax.set_xlim(self.t_start, self.t_stop) 
            ax2.set_ylim(0, 360) 
            ax.set_xlabel('Time [' + self.timeUnit + ']')
            ax.set_ylabel('Signal')
            ax2.set_ylabel('Phase [deg]')
            
            pl.show()
        
        if justPhase:
            return self.hilbertPhase
        else:
            return self.hilbertPhase, self.hilbertAbsolute, self.analyticalTrain

    
    def correlate(self, id0, id1, sigma=2., display=False):
        """ Cross-correlation of two spike trains.
        
        Parameters:
        id0, id1 ... IDs of the spiketrains
        display ... display the cross-correlogram
        
        Returns:
        tau ... time shift of the correlation
        corr ... correlation value
        
        """
        dt = .1
        xmax = 4.*sigma
        xmin = -xmax
        offset = -xmin / dt
        kernel = numpy.exp(-numpy.arange(xmin, xmax, dt) ** 2 / (2 * sigma ** 2))
        self.convolveSpikes(kernel, offset, dt, display=True)
        
        train0 = self.convolvedSpikeTrains[id0]
        train1 = self.convolvedSpikeTrains[id1]
        corr = numpy.correlate(train1, train0, 'same')  # correlation with respect to train0
        # tau = numpy.arange(-corr.shape[0]/2+1, corr.shape[0]/2+1)*dt
        tau = numpy.arange(-corr.shape[0] / 2, corr.shape[0] / 2) * dt
        if display:
            fig = pl.figure(103, figsize=(12, 7))
            ax = fig.add_subplot(111)   
            ax.plot(tau, corr, '-')
            ax.set_xlabel('Offset [' + self.timeUnit + ']')
            pl.show()
        
        return tau, corr

    def getSpikePhases(self):
        """
        Calculate the spike phases from the hilbert phase and attach this information
        to the particular spiketrains in the spikezugList.
        """
        
        try: self.hilbertPhase
        except AttributeError:
            self.hilbertTransform()
        
        for num in self.spiketrains:
            indices = numpy.array([(numpy.abs(self.timeAxis - s)).argmin() for s in self[num].spike_times])
            for index in numpy.where(indices >= self.timeAxis.shape[0] - 1)[0]:
                indices[index] = self.timeAxis.shape[0] - 1
            self[num].spike_phases = self.hilbertPhase[indices]


    def spectrogram(self, windowSize=4096, overlap=None, display=False, minFreq=None, maxFreq=None):
            
            try: self.filteredSummedConvolvedSpikeTrains
            except AttributeError:
                self.filterSummedConvolvedSpikeTrains()
            
            
            # parameters of the spectrogram
            if not overlap:
                overlap = int(windowSize * .9)
            samplingFreq = 1. / self.binwidth
            
            # use matplotlibs specgram function for the spectrogram
            Pxx, freqs, t = matplotlib.mlab.specgram(self.summedConvolvedSpikeTrains, \
                NFFT=windowSize, noverlap=overlap, pad_to=windowSize * 10, Fs=samplingFreq)           
            Pxxfiltered, freqs_filtered, t_filtered = matplotlib.mlab.specgram(self.filteredSummedConvolvedSpikeTrains, \
                NFFT=windowSize, noverlap=overlap, pad_to=windowSize * 10, Fs=samplingFreq)
            
            freqs *= 1000
            freqs_filtered *= 1000
            
            if display:
                if not minFreq:
                    minFreq = freqs.min()
                if not maxFreq:
                    maxFreq = freqs.max()
                
                indexstart = numpy.where(freqs >= minFreq)[0][0]
                indexend = numpy.where(freqs <= maxFreq)[0][-1]
                print indexstart, indexend
                
                fig = pl.figure(102, figsize=(10, 7))
                ax = fig.add_subplot(121)
                t, freqs = numpy.meshgrid(t, freqs)
                ax.pcolormesh(t[indexstart:indexend], freqs[indexstart:indexend], Pxx[indexstart:indexend])
                ax.set_ylim(minFreq, maxFreq)
                ax.set_xlabel('Time [' + self.timeUnit + ']')
                ax.set_ylabel('Frequency [Hz]')
                
                ax = pl.subplot(122)
                t_filtered, freqs_filtered = numpy.meshgrid(t_filtered, freqs_filtered)
                ax.pcolormesh(t_filtered[indexstart:indexend], freqs[indexstart:indexend], Pxxfiltered[indexstart:indexend])
                ax.set_ylim(minFreq, maxFreq)
                ax.set_xlabel('Time [ms]')
                
                pl.show()
            
            return (Pxx, freqs, t), (Pxxfiltered, freqs_filtered, t_filtered)

    
    def time_offset(self, offset):
        super(spikezugList, self).time_offset(offset)
        if hasattr(self, 'timeAxis'):
            self.timeAxis += offset

    def time_slice(self, t_start, t_stop, meta=None):
        """ Return a new spikezugList obtained by slicing between t_start and t_stop
        
        Parameters:
            t_start - begining of the new SpikeTrain, in ms.
            t_stop  - end of the new SpikeTrain, in ms.
        """
        if not meta:
            meta = self.meta
        
        dims = 1
        if meta.has_key('dimensions'):
           dims = meta['dimensions'] 
        new_spikezugList = spikezugList(t_start=t_start, t_stop=t_stop, dims=dims)
        new_spikezugList.meta = meta
        if inspect.ismethod(self.id_list):
            IDlist = self.id_list()
        else:
            IDlist = self.id_list
            
        for id in IDlist:
            new_spikezugList.append(id, self.spiketrains[id].time_slice(t_start, t_stop))
            new_spikezugList._SpikeList__calc_startstop()
            
        return new_spikezugList
    
