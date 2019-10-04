"""
signale.cscs
============

A module for signals from continously sampled channels (CSCs).
"""
__author__ = ("KT")
__version__ = "4.11, October 2013"

# python modules

# other modules
import numpy
import pylab
import scipy.signal as scsig
import NeuroTools.signals as NTsig
import matplotlib.pyplot as pl
import matplotlib as mpl
import spectrum as sp
from scipy import stats
import statsmodels.api as sm
import spectrum as sp
# custom made modules
import custom_plot

# package modules
from tools import findNearest
import signals
from signals import signal
# #import io
from sklearn import linear_model

###################################################### FUNCTIONS

###################################################### CLASSES

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

        factor = signals._getUnitFactor(self.timeUnit, newUnit)

        # change the times
        self.t_start *= factor
        self.t_stop *= factor
        self.dt *= factor

        self.timeUnit = newUnit


    def fft(self):
        '''
        function to calculate Fourier transform ussing fast Fourier transform implemented in pylab.
        :param display:
        :return:
        '''
        n = self.signal.shape[0]
        if self.timeUnit == 'ms':
            binwidth = self.dt / 1000.
        elif self.timeUnit == 's':
            binwidth = self.dt
        self.sp = pylab.rfft(self.signal)
        self.sp /= self.signal.size

    def periodogram(self,method='fft',display=False):
        nfft = sp.nextpow2(self.signal.size)
        if method=='fft':
            p = sp.Periodogram(self.signal,sampling=self.sampleFreq,NFFT=2**nfft)
            p();
        self.freq = numpy.array(p.frequencies())
        self.spPower = p.psd
        if display:
            p.plot()

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
        pl.yscale('log')
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
        self.cut_freqz = [minFreq,maxFreq]

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


    def spectrogram(self, minFreq=None, maxFreq=None, windowSize=2*4096, overlap=None, display=False,whiten=False):
        """
        Desciption is missing!!!
        """

        timeUnit = self.timeUnit
        if self.timeUnit != 'ms':
            self.changeTimeUnit('ms')

        # parameters of the spectrogram
        if not overlap:
            overlap = int(windowSize * .9)
        samplingFreq = 1. / self.dt

        # use matplotlibs specgram function for the spectrogram
        if whiten:
            if hasattr(self,'signal_white'):
                Pxx, freqs, t = mpl.mlab.specgram(self.signal_white, \
                                                  NFFT=windowSize, noverlap=overlap, pad_to=windowSize * 3, Fs=samplingFreq)
            else:
                print 'Signal is not whitened yet, whitening using ARMA(2,0) model by default...'
                self.whitenARMA()
                Pxx, freqs, t = mpl.mlab.specgram(self.signal_white, \
                                                  NFFT=windowSize, noverlap=overlap, pad_to=windowSize * 3, Fs=samplingFreq)
        else:
            Pxx, freqs, t = mpl.mlab.specgram(self.signal, \
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
            ax.set_xlabel('Time (' + self.timeUnit + ')')
            ax.set_ylabel('Frequency (Hz)')

            pl.show()

        # switch back time unit if necessary
        if self.timeUnit != timeUnit:
            self.changeTimeUnit(timeUnit)

        return Pxx, freqs, t


    def time_slice(self, t_start, t_stop):
        """ 
        Slice NeuralynxCSC between t_start and t_stop
        Parameters
        ----------
            t_start - begining of the new NeuralynxCSCList, in ms.
            t_stop  - end of the new NeuralynxCSCList, in ms.
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
            
        """

        assert self.t_start <= t_start
        assert self.t_stop >= t_stop

        t = self.time_axis()
        index1 = findNearest(t, t_start)[0]
        index2 = findNearest(t, t_stop)[0]
        self.signal = self.signal[index1:index2]
        self.recalc_timeAxis()
        
    def rms(self, window_size):
        '''
        Calculates the rms of the signal and the filtered signal!!!
        added by ACh Oct 2013
        
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        ssq = numpy.power(self.signal, 2)
        window = numpy.ones(window_size) / float(window_size)
        self.rms_signal = numpy.sqrt(numpy.convolve(ssq, window, 'same'))
        if hasattr(self, 'signal_filtered'):
            ssqf = numpy.power(self.signal_filtered, 2)
            window = numpy.ones(window_size) / float(window_size)
            self.rms_signal_filtered = numpy.sqrt(numpy.convolve(ssqf, window, 'same'))
        else:
            print 'Currently there is no filtered signal. \n \
            if you have filtered the signal and you want the RMS of filtered \n \
            signal you need to call this function once more after filteration!'
    def whitenARMA(self,AR=2,MA=0):
        '''
        whiten the signal using ARMA!
        :param AR:
        :param MA:
        :return:
        '''
        arma = sm.tsa.ARMA(self.signal, (AR,MA)).fit(disp=0)
        print 'ARMA parameters calculated for order(%s,%s)' %(AR,MA)
        self.signal_white = arma.resid
        self.AR_order = AR
        self.MA_order = MA


    def whitenLinearRegression(self):
        '''
        whiten the signal using flat psd obtained by estimateNoiseFreqRelatyion
        :return:
        '''
        if not hasattr(self,'spPower_whitened'):
            self.flattenPowerspectrum()
        self.signal_whiteLS = numpy.fft.ifft(self.sp)
        self.signal_whiteLS *= self.signal.size  # rescale from normalization of the fourier components
        self.signal_whiteLS = self.signal_whiteLS[:self.signal.size]


    def estimateNoisefreqRelation(self,detectOutliers = True,noOfIterations = 0,display=False):
        '''
        a function to estimate the alpha in Sp(N)~1/f^alpha
        Parameters:
        detectoutliers: if True, the exponent is estimated using RANSAC method.
        Otherwise normal linear regression without considering outliers is used.

        :return:
         noiseFreqExponent:  the power,frequency relation exponent!

         References:
         https://en.wikipedia.org/wiki/Colors_of_noise
         http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.RANSACRegressor.html
        '''
        if not hasattr(self,'spPower'):
            self.fft()
        cnd = self.freq>10
        ff = numpy.log10(numpy.array([self.freq[cnd]]).T)
        pp = numpy.log10(self.spPower[cnd])
        model = linear_model.LinearRegression()
        model.fit(ff,pp)
        # Robustly fit linear model with RANSAC algorithm
        model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
        model_ransac.fit(ff,pp)
        inlier_mask = model_ransac.inlier_mask_
        outlier_mask = numpy.logical_not(inlier_mask)
        if display:
            # Predict data of estimated models
            line_X = numpy.arange(.5, 3.54)
            line_y = model.predict(line_X[:, numpy.newaxis])
            line_y_ransac = model_ransac.predict(line_X[:, numpy.newaxis])
            pl.plot(ff[inlier_mask], pp[inlier_mask], '.g', label='Inliers')
            pl.plot(ff[outlier_mask], pp[outlier_mask], '.r', label='Outliers')
            pl.plot(line_X, line_y, '-k', label='Linear regressor')
            pl.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
            pl.legend(loc='lower right')
            pl.show()
        if detectOutliers:
            if noOfIterations:
                for ii in range(noOfIterations):
                    coefList = []
                    model_ransac.fit(ff,pp)
                    coefList.append(model_ransac.estimator_.coef_)
                    self.noiseFreqExponent = numpy.mean(coefList)
                else:
                    self.noiseFreqExponent = model_ransac.estimator_.coef_[0][0]
        else:
            self.noiseFreqExponent = model.coef_[0]
    def flattenPowerspectrum(self):
        '''
        A function to whiten the power spectrum using the estimated noise color.
        Parameters
        ----------


        Returns
        ----------

        See also
        ----------
         estimateNoisefreqRelation
        Notes
        ----------
        https://en.wikipedia.org/wiki/Colors_of_noise
        '''
        if not hasattr(self,'noiseFreqExponent'):
            self.estimateNoisefreqRelation()

        freqz = numpy.power(numpy.array(self.freq),-1.0*self.noiseFreqExponent)
        self.spPower_whitened = numpy.multiply(self.spPower,freqz)


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
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """

        factor = signals._getUnitFactor(self.timeUnit, newUnit)

        for csc in self:
            csc.changeTimeUnit(newUnit)

        self.t_start *= factor
        self.t_stop *= factor
        self.dt *= factor
        self.timeUnit = newUnit

    def getLimits(self):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        
        self.min = 1e100
        self.max = -1e-100
        for csc in self:
            self.min = min(self.min, csc.signal.min())
            self.max = max(self.max, csc.signal.max())

    def plot(self, fig=None,ax=None):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        try: self.min
        except AttributeError:
            self.getLimits()
        if not fig:
            fig = pl.figure(figsize=(12, 7))
        if not ax:
            ax = fig.add_subplot(111)
        ylabels = []
        last_id = self.id_list()[-1] + 1
        for id, csc in enumerate(self):
            s = csc.signal.copy()
            s -= s.mean()
            s /= self.max * .8
            ax.plot(csc.time_axis(), s + id, '-', linewidth=1)
            ylabels.append(csc.tags['channel'].split('.')[0])
        custom_plot.huebschMachen(ax)
        ax.set_xlabel('Time (' + self.timeUnit + ')')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, rotation=0)
        ax.set_xlim(self.t_start, self.t_stop)
        ax.set_ylim(-1, id + 1)

        return fig, ax


    def fft_plot(self, freq_min=0., freq_max=60.):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
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
        ax.set_xlabel('Frequency (Hz)')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, rotation=0)
        ax.set_xlim(freq_min, freq_max)
        ax.set_ylim(0, id + 1)

        return fig, ax

    def removeMean(self):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        for id in self.analog_signals:
            self[id].removeMean()

    def showTags(self):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        signal.showTags(self)
        for id in self.analog_signals:
            print id
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
        self.ADBitVolts = 1
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
                self.timeAxis = numpy.append(self.timeAxis, numpy.linspace(t, self.times[i + 1], self.tags['validSamples'][i], endpoint=False))
            self.timeAxis = numpy.append(self.timeAxis, numpy.arange(self.times[-1], self.times[-1] + self.dt * (self.tags['validSamples'][-1]), self.dt))
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
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        factor = signals._getUnitFactor(self.timeUnit, newUnit)

        # change the times
        self.t_start *= factor
        self.t_stop *= factor
        self.dt *= factor

        self.times_start *= factor
        self.times_stop *= factor
        self.times *= factor
        self.timeAxis *= factor

        self.timeUnit = newUnit
    def theta_delta_ratio(self, time_window=0.5,whiten =True ):
        '''
        This function calculates the spectrograms in theta and delta range and returns the ratio of them.
        
        Parameters:
                time_window ... is the time window (in seconds) over which the spectrogram is calculated.
        '''
        window_size = int(time_window / self.dt * 1e3)
        theta, freq, t = self.spectrogram(minFreq=5, maxFreq=10, windowSize=window_size,whiten=whiten)
        delta, freq, t = self.spectrogram(minFreq=0.05, maxFreq=4.9, windowSize=window_size,whiten=whiten)
        self.th_del_ratio = [t[1, :], (numpy.sum(theta, axis=0) / numpy.sum(delta, axis=0))]
        
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

        self.analyticalTrain = scsig.hilbert(self.signal_filtered.real)
        self.hilbertAbsolute = numpy.absolute(self.analyticalTrain)
        self.hilbertPhase = numpy.angle(self.analyticalTrain, deg=True)
        self.hilbertAbsSmooth = numpy.convolve(self.hilbertAbsolute, scsig.gaussian(100, 33), 'same')
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
            ax.set_xlabel('Time (' + self.timeUnit + ')')
            ax.set_ylabel('Absolute value')
            ax2.set_ylabel('Phase (deg)')

            ax = fig.add_subplot(313)
            ax2 = ax.twinx()
            ax.plot(self.timeAxis, self.analyticalTrain.real)
            ax2.plot(self.timeAxis, self.hilbertPhase, 'g')
            ax.set_xlim(self.t_start, self.t_stop)
            ax2.set_ylim(0, 360)
            ax.set_xlabel('Time (' + self.timeUnit + ')')
            ax.set_ylabel('Signal')
            ax2.set_ylabel('Phase (deg)')

            pl.show()

        if justPhase:
            return self.hilbertPhase
        else:
            return self.hilbertPhase, self.hilbertAbsolute, self.analyticalTrain


    def plot(self, fig=None, ax=None):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        if not fig:
            fig = pl.figure(figsize=(12, 7))
        if not ax:
            ax = fig.add_subplot(111)

        ax.plot(self.time_axis(), self.signal, '-')

        ax.set_xlabel('Time (' + self.timeUnit + ')')
        ax.set_xlim(self.t_start, self.t_stop)

        pl.show()

        return fig, ax


    def purge(self):
        """ Remove values that fill up the possible range.
        
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        
        """

        ds = numpy.diff(self.signal)
        indices, = numpy.where(ds == 0)
        indices = numpy.concatenate((indices - 1, indices, indices + 1))
        self.signal = numpy.delete(self.signal, indices)

        # to do properly:
        # self.times ... purge as well

        self.recalc_timeAxis()

    def resample(self, new_dt):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
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
        
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
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

        index1 = numpy.searchsorted(self.timeAxis, t_start)
        index2 = numpy.searchsorted(self.timeAxis, t_stop)
        print index1,index2
        if index2 > self.timeAxis.shape[0]:
            index2 = self.times.shape[0]
        tAxis = self.timeAxis[index1:index2]
        slicedSignal = self.signal[index1:index2]
        #print tAxis,slicedSignal
        slicedChunk = NeuralynxCSC(tAxis,slicedSignal,sampleFreq=self.sampleFreq)
        if hasattr(self,'signal_filtered'):
            slicedChunk.signal_filtered = self.signal_filtered[index1:index2]
        if hasattr(self,'signal_white'):
            slicedChunk.signal_white = self.signal_white[index1:index2]
        if hasattr(self,'signal_purged'):
            slicedChunk.signal_purged = self.signal_purged[index1:index2]
        self.recalc_timeAxis()
        return slicedChunk
    def REM_detector(self, cutThresh=None,peakTresh=None):
        '''
        This function detects the REM Sleep episodes.
        
        Parameters
        ----------
            thresh ... Threshold on theta/delta ratio which is basically the selection criteria for REM episodes. The default value is (Mean + 1sigma)
             of theta/delta ratio.
        
        Returns
        ----------
            rem_episodes... A numpy array of shape(n,2) containing the t_start and t_stop of each REM episode on each row.
            sws_episodes... A numpy array of shape(n,2) containing the t_start and t_stop of each SWS episode on each row.
        
        See also
        ----------
        REM_signal, theta_delta_ratio, 
        Notes
        ----------
        
        '''
        if hasattr(self,'th_del_ratio') is False:
            print 'theta/delta ratio is not calculated yet.\n Calculating th/del ratio...'
            self.theta_delta_ratio()
            
        ratio = self.th_del_ratio[1]
        avg = ratio.mean()
        std = ratio.std()
        if  cutThresh:
            cutThresh = cutThresh
        else:
            cutthresh = avg
        if peakTresh:
            peakTresh = peakTresh
        else:
            peakTresh = avg + std

        rem = numpy.where(ratio > cutthresh)[0]
        rem_break = numpy.where(numpy.diff(rem) > 1)[0]
        rem_split = numpy.split(rem, rem_break + 1)
        rem_times = []
        sws_times = []
        for item in rem_split:
            if item.size > 1 and ratio[item].max() > peakTresh:
                rem_times.append([self.th_del_ratio[0][item[0]], self.th_del_ratio[0][item[-1]]])
        rem_times = numpy.array(rem_times)
        first = numpy.array([])
        second = numpy.array([])
        first = numpy.append(first, [item[1] for item in rem_times])
        second = numpy.append(second, [item[0] for item in rem_times])
        first = numpy.insert(first, 0, 0)
        second = numpy.append(second, self.time_axis()[-1])
        for ii in range(first.size):
            sws_times.append([first[ii], second[ii]])
        sws_times = numpy.array(sws_times)
        self.rem_episodes = rem_times
        self.sws_episodes = sws_times
        
    def REM_signal(self):
        ''' 
        This function selects the REM episodes from a csc signal!
        
        Parameters
        ----------
        
        
        Returns
        ----------
        rem_times:
        rem_signal:
        See also
        ----------
        REM_detector, theta_delta_ratio,
        Notes
        ----------
        '''
        rem_time = []
        rem_signal = []
        t0 = self.time_axis()
        for item in self.rem_episodes:
            rem_idx = numpy.where((t0 > item[0]) & (t0 < item[1]) == True)[0]
            rem_time.append(t0[rem_idx])
            rem_signal.append(self.signal[rem_idx])
        self.rem_times = numpy.array(rem_time)
        self.rem_signal = numpy.array(rem_signal)

    def SWS_signal(self,filter = False,f_1 = None,f_2 = None):
        '''
        This function selects the SWS episodes from a csc signal!
        
        Parameters
        ----------
        filter:
        f_1:
        f_2:
        
        Returns
        ----------
        sws_signal:
        sws_times:
        
        See also
        ----------
        REM_detector, theta_delta_ratio,
        Notes
        ----------
        
        '''
        if hasattr(self,'sws_episodes')is False:
            print 'SWS/REM episodes are not detected yet!\n Detecting using REM_Detector function with default thershold value... '
            self.REM_detector()
        sws_time = []
        sws_signal = []
        sws_signal_filtered = []
        t0 = self.timeAxis
        for item in self.sws_episodes:
            sws_idx = numpy.where((t0 > item[0]) & (t0 < item[1]) == True)[0]
            sws_time.append(t0[sws_idx])
            sws_signal.append(self.signal[sws_idx])
            if filter:
                sws_signal_filtered.append(self.signal_filtered[sws_idx])
        self.sws_times = numpy.array(sws_time)
        self.sws_signal = numpy.array(sws_signal)
        if filter:
            self.sws_filtered = numpy.array(sws_signal_filtered)
    
    
    def ripple_recorder(self, sigma=20, length=20, rippleMix=False, rippleCut=True,SWRmix=True,removeREMripples = False):
        '''
            Function for detecting and recording ripples in a single csc! It takes a csc object,
            detect the ripples and add them as an attribute to the csc object.
        
        Parameters
        ----------
            sigma... sigma  of gaussian of window of length (lenght*sigma) to be used for smoothing the hilbertAbs!
            length... length of gausssian window in units of sigma!
            rippleMix...  if True, function will glue all the detected ripples and return them as rippMX attribute!
            rippleCut...  if True, function will make an array with the  same size of signal , only containing ripples with
                          Everything else is set to zero!
        
        Returns
        ----------
            ripples
            rippMX
            rippCut
        See also
        ----------
        
        Notes
        ----------
            
                               
    
            Last Modified June 2014
        '''
        ripples = []
        if hasattr(self, 'signal_filtered') is False:
            print 'Signal is not filtered yet!\n Filtering the signal with the default values...'
            self.filter(100, 250)
        if hasattr(self, 'hilbertAbsolute') is False:
            print 'Hilbert transform is not calculated!\n Calculating Hilbert transform of the signal...'
            self.hilbertTransform()
        if hasattr(self,'sws_signal') is False:
            print 'SWS and REM signals are not trimmed yet!\n Trimming SWS/REM signals...'
            self.SWS_signal(filter=True,f_1=100,f_2=250)
        sws_hilbert_smooth =[]       #This will keep the hilbertAbsolute smoothed signal in separate chunks!!!
        for ii,item in enumerate(self.sws_filtered):
            print '%0.3f prcnt ---' %(1e2*ii/len(self.sws_filtered))
            if item.size:
                #print item.size
                sws_hilbert_smooth.append(numpy.convolve(numpy.absolute(scsig.hilbert(item)),scsig.gaussian(length * sigma,sigma), 'same'))
        self.sws_hilbert_smooth = sws_hilbert_smooth
        print 'Sir Hilbert is Abs and Smooooooooooooth!!!'
        sws_hilbert_smooth = numpy.hstack(sws_hilbert_smooth)
        print '#'
        smooth = numpy.convolve(self.hilbertAbsolute, scsig.gaussian(length * sigma, sigma), 'same')
        print '###'
        std = sws_hilbert_smooth.std()
        avg = sws_hilbert_smooth.mean()
        print '#####'
        sig1 = numpy.where(smooth > (avg + std))[0]
        sig1_break = numpy.where(numpy.diff(sig1) > 1)[0]
        sig1_split = numpy.split(sig1, sig1_break + 1)
        print '#######'
        t_chain = []
        t_start = 0.0
        t_end = 0.0
        peak = 0.0
        t_peak = 0.0
        ripple_mix = []
        SWR_mix = []
        ripple_cut = numpy.zeros(self.signal.size)
        print sig1_split
        for i in range(len(sig1_split)):
            print sig1_split[i]
            if  smooth[sig1_split[i]].max() > avg + 3 * std:
                t_chain = numpy.where(smooth == smooth[sig1_split[i]].max())[0]
                t_start = self.timeAxis[sig1_split[i][0]]
                t_end = self.timeAxis[sig1_split[i][-1]]
                peak = smooth[sig1_split[i]].max()
                t_peak = self.timeAxis[numpy.intersect1d(sig1_split[i], t_chain)][0]
                ripples.append([t_start, t_end, peak, t_peak])
            if rippleMix == True:
                ripple_mix.append(self.signal_filtered[sig1_split[i]])
            if SWRmix == True:
                SWR_mix.append(self.signal[sig1_split[i]])
            if rippleCut == True:
                ripple_cut[sig1_split[i]] = self.signal_filtered[sig1_split[i]]   
        if rippleMix == True:
            ripple_mix = numpy.array(ripple_mix)
            self.rippMX = ripple_mix
        if rippleCut == True:
            self.signal_just_ripples = ripple_cut
        if SWRmix == True:
            SWR_mix = numpy.array(SWR_mix)
            self.SWRmix = SWR_mix
        ripples = numpy.array(ripples)
        self.ripples = ripples
        self.hilbertAbsSmooth = smooth
        print 'Ripple detection on' , self.tags['file'], ' is DONE!!!\n \
        Now you should see the attribute ripples with [t_start, t_end,peak value, t_peak] in each row for detected SWRs.'
        if removeREMripples:
            SUM = 0
            j = 0
            for item in self.rem_episodes:
                for jtem in self.ripples:
                    peak_in = jtem[-1] <item[1] and jtem[-1]>item[0]
                    starts_in = jtem[0] <item[1] and jtem[0]>item[0]
                    ends_in = jtem[1] <item[1] and jtem[1]>item[0]
                    if peak_in or starts_in or ends_in :
                        ripples = numpy.delete(ripples,j-SUM,0)
                        SUM += 1
                        j+= 1
        self.SWS_ripples = ripples
                    
    def ripple_statistics(self,fig = None, ax1= None,ax2=None,ax3= None,ax4= None, sws = True):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        if not hasattr(self, 'ripples'):
            print "Ripples are not detected yet! "
        else:
            duration = []
            tresh =[]
            if sws:
                for item in self.SWS_ripples:
                    duration.append(item[1] - item[0])
            else:
                for item in self.ripples:
                    duration.append(item[1] - item[0])
            
        for item in self.sws_hilbert_smooth:
            tresh.append(item.mean()+3*item.std())
        tresh = numpy.array(tresh)
        if not (fig or ax1 or ax2 or ax3 or ax4):
            fig = pl.figure()
            ax1 = fig.add_subplot(221)
            ax2 = fig.add_subplot(222)
            ax3 = fig.add_subplot(223)
            ax4 = fig.add_subplot(224)
        #####################################################
        if ax1:
            ax1.hist(duration, bins=100)
            ax1.set_title('ripp duration')
            ax1.locator_params(nbins=4)
            ax1.set_xlim(20,150)
        ####################################################
        if ax2: 
            ax2.scatter(range(tresh.size),tresh,marker=(5,0))
            ax2.axhline(tresh.mean(),linewidth=1, color='r', alpha=0.7, label='Mean',
                    linestyle='--')
            ax2.axhline(self.hilbertAbsSmooth.mean()+3*self.hilbertAbsSmooth.std(),linewidth=1, color='g', alpha=0.7, label='Mean,Signal',
                    linestyle='-.')
            ax2.axhline(numpy.hstack(self.sws_hilbert_smooth).mean()+3*numpy.hstack(self.sws_hilbert_smooth).std(),linewidth=1, color='k', alpha=0.7, label='Mean,Signal',
                    linestyle='--')
            ax2.set_title('Thresholds')
            ax2.locator_params(nbins=4)
        ####################################################
        if ax3:
            if not hasattr(self, 'IRI'):
                self.iri()
            ax3.hist(self.IRI, bins=400)
            ax3.set_title('IRI')
            ax3.locator_params(nbins=4)
            ax3.set_xlim(0,4000)
        ####################################################
        #doing the FFT!
        if ax4:
            self.ripple_fft(3)
            fft_smth = numpy.convolve(self.ripp_fft_raw[1],scsig.gaussian(700,200),'same')
            ax4.fill_between(self.ripp_fft_raw[0],fft_smth,facecolor='yellow',alpha=0.4)
            ax4.plot(self.ripp_fft_raw[0],fft_smth,alpha=0.7)
            idx = (numpy.abs(self.ripp_fft_raw[0]-150)).argmin()
            ax4.set_xlim(100,220)
            ax4.set_ylim(0,2 * fft_smth[idx])
        pl.tight_layout()
        return fig,ax1,ax2,ax3,ax4

    def rmsPlot(self):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        mv_converter = self.ADBitVolts * 1e6  # converting factor to mV for the ampilitudes
        f, axarr = pl.subplots(4, sharex=True)
        # #Hilbert transformed signal
        axarr[0].plot(self.timeAxis[0:self.signal.size] / 1000, self.hilbertAbsolute * mv_converter)
        axarr[0].axhline(y=(self.hilbertAbsolute.mean() + self.hilbertAbsolute.std()) * mv_converter, color='k')
        axarr[0].axhline(y=(self.hilbertAbsolute.mean() + 3 * self.hilbertAbsolute.std()) * mv_converter, color='g')
        axarr[0].set_title('Hilbert Detection Method')
        axarr[0].set_xlim(self.timeAxis[0] / 1000, self.timeAxis[-2] / 1000)
        
        # #filtered_rms signal 
        axarr[1].plot(self.timeAxis[0:self.signal.size] / 1000, self.rms_signal_filtered * mv_converter, 'm')
        axarr[1].axhline(y=(self.rms_signal_filtered.mean() + self.rms_signal_filtered.std()) * mv_converter, color='k')
        axarr[1].axhline(y=(self.rms_signal_filtered.mean() + 3 * self.rms_signal_filtered.std()) * mv_converter, color='g')
        axarr[1].set_title('RMS Detection Method')
        
        # #filtered signal
        axarr[2].plot(self.timeAxis[0:self.signal.size] / 1000, self.signal_filtered * mv_converter, 'r')
        axarr[2].set_ylabel('Ampilitude(mV)')
        axarr[2].set_title('Filtered')
        
        # #raw signal
        axarr[3].plot(self.timeAxis[0:self.signal.size] / 1000, self.signal * mv_converter, 'k')
        axarr[3].set_xlabel('Time(s)')
        axarr[3].set_title('Signal')
        f.suptitle(self.tags.get('file') + 'on channel: ' + str(self.tags.get('channel')))
        pl.show()

    def ripplePlot(self,sws_ripples=True):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        f, axarr = pl.subplots(4, sharex=True,figsize=[20,15])
        if not hasattr(self, 'ripples'):
            print 'Ripples are not detected yet, please use ripplerecorder function to detect them first!'
        elif sws_ripples:
            mv_converter = self.ADBitVolts * 1e6  # converting factor to mV for the ampilitudes
            # #filtered_HilbertAbsolute signal 
            for i in range(len(self.SWS_ripples)):
                axarr[0].axvspan(self.SWS_ripples[i][0], self.SWS_ripples[i][1], facecolor='g', alpha=0.5)
                axarr[0].axvline(self.SWS_ripples[i][-1], linewidth=1)
            axarr[0].plot(self.timeAxis[0:self.signal.size], self.hilbertAbsolute * mv_converter, 'm')
            axarr[0].set_title('HilbertAbs of filtered signal')
            axarr[0].axhline(y=(self.hilbertAbsolute.mean() + self.hilbertAbsolute.std()) * mv_converter, color='k', linewidth=1)
            axarr[0].axhline(y=(self.hilbertAbsolute.mean() + 3 * self.hilbertAbsolute.std()) * mv_converter, color='g', linewidth=1)
            axarr[0].set_xlim(self.timeAxis[0], self.timeAxis[-2])
            axarr[0].set_yticks([]) 
            
            # #filtered_HilbertAbsolute signal 
            for i in range(len(self.SWS_ripples)):
                axarr[1].axvspan(self.SWS_ripples[i][0], self.SWS_ripples[i][1], facecolor='g', alpha=0.5)
                axarr[1].axvline(self.SWS_ripples[i][-1], linewidth=1)
            axarr[1].plot(self.timeAxis[0:self.signal.size], self.hilbertAbsSmooth,)
            axarr[1].set_title('HilbertAbsSmooth of filtered signal')
            axarr[1].axhline(y=(self.hilbertAbsSmooth.mean() + self.hilbertAbsSmooth.std()), color='k', linewidth=1)
            axarr[1].axhline(y=(self.hilbertAbsSmooth.mean() + 3 * self.hilbertAbsSmooth.std()), color='g', linewidth=1)
            axarr[1].set_xlim(self.timeAxis[0], self.timeAxis[-2])
            axarr[1].set_yticks([])
            
            # #filtered signal
            for i in range(len(self.SWS_ripples)):
                axarr[2].axvspan(self.SWS_ripples[i][0], self.SWS_ripples[i][1], facecolor='g', alpha=0.5)
            axarr[2].plot(self.timeAxis[0:self.signal.size], self.signal_filtered * mv_converter, 'r')       
            axarr[2].set_yticks([])
            axarr[2].set_title('Filtered Signal')
            
                        
            # #raw signal
            for i in range(len(self.SWS_ripples)):
                axarr[3].axvspan(self.SWS_ripples[i][0], self.SWS_ripples[i][1], facecolor='g', alpha=0.5)
            axarr[3].plot(self.timeAxis[0:self.signal.size], self.signal * mv_converter, 'k',lw=0.5)
            axarr[3].set_xlabel('Time(ms)')
            axarr[3].set_title('Raw Signal')
            axarr[3].set_yticks([])
            f.suptitle(self.tags.get('file') + ' on channel: ' + str(self.tags.get('channel')))
            pl.show()
        else:
            mv_converter = self.ADBitVolts * 1e6  # converting factor to mV for the ampilitudes
            # #filtered_HilbertAbsolute signal
            for i in range(len(self.ripples)):
                axarr[0].axvspan(self.ripples[i][0], self.ripples[i][1], facecolor='g', alpha=0.5)
                axarr[0].axvline(self.ripples[i][-1], linewidth=1)
            axarr[0].plot(self.timeAxis[0:self.signal.size], self.hilbertAbsolute * mv_converter, 'm')
            axarr[0].set_title('HilbertAbs of filtered signal')
            axarr[0].axhline(y=(self.hilbertAbsolute.mean() + self.hilbertAbsolute.std()) * mv_converter, color='k', linewidth=1)
            axarr[0].axhline(y=(self.hilbertAbsolute.mean() + 3 * self.hilbertAbsolute.std()) * mv_converter, color='g', linewidth=1)
            axarr[0].set_xlim(self.timeAxis[0], self.timeAxis[-2])
            axarr[0].set_yticks([])

            # #filtered_HilbertAbsolute signal
            for i in range(len(self.ripples)):
                axarr[1].axvspan(self.ripples[i][0], self.ripples[i][1], facecolor='g', alpha=0.5)
                axarr[1].axvline(self.ripples[i][-1], linewidth=1)
            axarr[1].plot(self.timeAxis[0:self.signal.size], self.hilbertAbsSmooth,)
            axarr[1].set_title('HilbertAbsSmooth of filtered signal')
            axarr[1].axhline(y=(self.hilbertAbsSmooth.mean() + self.hilbertAbsSmooth.std()), color='k', linewidth=1)
            axarr[1].axhline(y=(self.hilbertAbsSmooth.mean() + 3 * self.hilbertAbsSmooth.std()), color='g', linewidth=1)
            axarr[1].set_xlim(self.timeAxis[0], self.timeAxis[-2])
            axarr[1].set_yticks([])

            # #filtered signal
            for i in range(len(self.ripples)):
                axarr[2].axvspan(self.ripples[i][0], self.ripples[i][1], facecolor='g', alpha=0.5)
            axarr[2].plot(self.timeAxis[0:self.signal.size], self.signal_filtered * mv_converter, 'r')
            axarr[2].set_yticks([])
            axarr[2].set_title('Filtered Signal')


            # #raw signal
            for i in range(len(self.ripples)):
                axarr[3].axvspan(self.ripples[i][0], self.ripples[i][1], facecolor='g', alpha=0.5)
            axarr[3].plot(self.timeAxis[0:self.signal.size], self.signal * mv_converter, 'k')
            axarr[3].set_xlabel('Time(ms)')
            axarr[3].set_title('Raw Signal')
            axarr[3].set_yticks([])
            f.suptitle(self.tags.get('file') + ' on channel: ' + str(self.tags.get('channel')))
            pl.show()

    def iri(self,sws=True):
        ''' 
        this function calculates the inter-ripple intervals (peak to peak) 
        written by AChenani
        
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        if not hasattr(self, 'ripples'):
            print 'Ripples are not detected yet, please use ripple_recorder() function to detect them first!'
        elif self.ripples.shape[0] < 3:
            print 'There are no enough ripples to calculate IRI!!!'
        elif sws:     
            iri = numpy.diff(self.SWS_ripples[:, 3])
            self.IRI = iri
        else:
            iri = numpy.diff(self.ripples[:, 3])
            self.IRI = iri
        return 
    
    def ripple_fft(self, method=None, d_thresh=None):
        '''
        This function calculates the frequency content of ripples.
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        ripp_mx = []
        if method == 1: 
            n = self.signal_just_ripples.size
            binwidth = self.dt / 1000
            sp = numpy.fft.rfft(self.signal_just_ripples)
            spPower = 2 * sp * numpy.conj(sp)
            freq = numpy.fft.fftfreq(n, d=binwidth)[0:n / 2 + 1]
            self.ripp_fft_zero_pad = [freq, spPower]
        if method == 2:
            if not d_thresh == None:
                for item in self.rippMX:
                    if len(item) > 2 * d_thresh:
                        ripp_mx.append(item)
            else:
                ripp_mx = self.rippMX
            rippmx = numpy.array([val for subl in ripp_mx for val in subl])
            rippmx = rippmx.flatten()
            n = rippmx.size
            binwidth = self.dt / 1000
            sp = numpy.fft.rfft(rippmx)
            spPower = 2 * sp * numpy.conj(sp)
            freq = numpy.fft.fftfreq(n, d=binwidth)[0:n / 2 + 1]
            self.ripp_fft_glue = [freq, spPower]
        if method == 3:
            swrmx = numpy.hstack(self.SWRmix) 
            n = swrmx.size
            binwidth = self.dt / 1000
            sp = numpy.fft.rfft(swrmx)
            spPower = 2 * sp * numpy.conj(sp)
            freq = numpy.fft.fftfreq(n, d=binwidth)[0:n / 2 + 1]
            self.ripp_fft_raw = [freq, spPower]
            print "SMILE"
            

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


    def plot(self, fig= None,ax=None):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """

        try: self.min
        except AttributeError:
            self.getLimits()

        if not fig:
            fig = pl.figure(figsize=(12, 7))
        if not ax:
            ax = fig.add_subplot(111)


        ylabels = []
        last_id = self.id_list()[-1] + 1
        for id, csc in enumerate(self):
            s = csc.signal.copy()
            s -= s.mean()
            s /= self.max * .8
            ax.plot(csc.time_axis(), s + id, '-', linewidth=1)
            ylabels.append(csc.tags['file'].split('.')[0])
        custom_plot.huebschMachen(ax)
        ax.set_xlabel('Time (' + self.timeUnit + ')')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, rotation=0)
        ax.set_xlim(self.t_start, self.t_stop)
        ax.set_ylim(-1, id + 1)

        return fig, ax


    def fft_plot(self, freq_min=0., freq_max=60.,fig= None,ax=None):
        """
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        """
        try: self.min
        except AttributeError:
            self.getLimits()

        if not fig:
            fig = pl.figure(figsize=(12, 7))
        if not ax:
            ax = fig.add_subplot(111)

        ylabels = []
        last_id = self.id_list()[-1] + 1
        #ax = fig.add_subplot(111)
        for id, csc in enumerate(self):
            if not hasattr(csc, 'spPower'):
                csc.fft()
            s = csc.spPower
            s -= s.mean()
            s /= s.max() * 1.1
            ax.plot(csc.freq, s + id, '-', linewidth=1)
            ylabels.append(csc.tags['file'].split('.')[0])
        custom_plot.huebschMachen(ax)
        ax.set_xlabel('Frequency (Hz)')
        yticks = numpy.arange(last_id)
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, rotation=0)
        #ax.set_xlim(freq_min, freq_max)
        ax.set_ylim(0, id + 1)

        return fig, ax

    def cscCompare(self, csc_no=0):
        '''
        This function plots all the raw signals with SWR events detected on one of the signals indicated by csc_no variable.
        Last Modified Oct. 2013
        
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        if csc_no > len(self) or csc_no < 1:
            print 'The argument must be between 1 and the length of cscList object!!!'
        elif not hasattr(self[csc_no - 1], 'ripples'):
            print 'Ripples are not detected yet, please use ripplerecorder function to detect them first!\n Please try agin...' 
        else:
            f, axarr = pl.subplots(len(self) + 2, sharex=True)
            mv_converter = self[csc_no - 1].ADBitVolts * 1e6  # converting factor to mV for the ampilitudes
            # #filtered_rms signal 
            for i in range(len(self[csc_no - 1].ripples)):
                axarr[0].axvspan(self[csc_no - 1].ripples[i][0] / 1000., self[csc_no - 1].ripples[i][1] / 1000., facecolor='g', alpha=0.5)
                axarr[0].axvline(self[csc_no - 1].ripples[i][-1] / 1000., linewidth=1)
            axarr[0].plot(self[csc_no - 1].timeAxis[0:self[csc_no - 1].signal.size] / 1000., self[csc_no - 1].hilbertAbsolute * mv_converter, 'm')
            axarr[0].axhline(y=(self[csc_no - 1].hilbertAbsolute.mean() + self[csc_no - 1].hilbertAbsolute.std()) * mv_converter, color='k', linewidth=1)
            axarr[0].axhline(y=(self[csc_no - 1].hilbertAbsolute.mean() + 3 * self[csc_no - 1].hilbertAbsolute.std()) * mv_converter, color='g', linewidth=1)
            axarr[0].set_xlim(self[csc_no - 1].timeAxis[0] / 1000, self[csc_no - 1].timeAxis[-2] / 1000)
            axarr[0].set_title(self[csc_no - 1].tags.get('file') + 'on channel: ' + str(self[csc_no - 1].tags.get('channel')), fontsize=10)

            # #filtered signal
            for i in range(len(self[csc_no - 1].ripples)):
                axarr[1].axvspan(self[csc_no - 1].ripples[i][0] / 1000., self[csc_no - 1].ripples[i][1] / 1000., facecolor='g', alpha=0.5)
            axarr[1].plot(self[csc_no - 1].timeAxis[0:self[csc_no - 1].signal.size] / 1000, self[csc_no - 1].signal_filtered * mv_converter, 'r')       
            axarr[1].set_ylabel('Ampilitude(mV)')
            axarr[1].set_title('Filtered')
            axarr[1].set_title(self[csc_no - 1].tags.get('file') + 'on channel: ' + str(self[csc_no - 1].tags.get('channel')), fontsize=10)

                        
            # #raw signal
            for k in range(2, len(axarr)):
                for i in range(len(self[k - 2].ripples)):
                    axarr[k].axvspan(self[k - 2].ripples[i][0] / 1000., self[k - 2].ripples[i][1] / 1000., facecolor='g', alpha=0.5)
                axarr[k].plot(self[k - 2].timeAxis[0:self[k - 2].signal.size] / 1000, self[k - 2].signal * mv_converter, 'k')
                # axarr[k].set_xlabel('Time(s)')
                axarr[k].set_title(self[k - 2].tags.get('file') + 'on channel: ' + str(self[k - 2].tags.get('channel')), fontsize=10)
            # f.suptitle(self[1].tags.get('file') + 'on channel: ' + str(self.tags.get('channel')))
            pl.show()
            return axarr
    
    def avg_csc(self):
        '''
        This function averages over all csc objects in cscList and adds this average as a new csc object to the end
        of cscList with proper tags!
        Last Modified Nov. 2013
        
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        avg = numpy.zeros(self[0].signal.size)
        for i in range(len(self)):
            avg = numpy.add(avg, self[i].signal)
        avg /= float(len(self))
        csc_avg = NeuralynxCSC(self[0].times, avg, self[0].sampleFreq)
        csc_avg.tags['file'] = 'CSC Average'
        csc_avg.tags['channel'] = 'ZDF'
        self.append(len(self), csc_avg)
        
    def rippleDetect(self, sigma=10, length=10):
        '''
            Function for detecting and recording ripples in an experiment! It takes a cscList object,
            averages the smoothed hilbert signal and then detect the ripples using threshhold method
            and add them as an attribute to the cscList object.
            Parameters
            ----------
            sigma: window length of the gaussian kernel to be used for smoothing the HilbertAbs signal.
            default value is zero i.e. no smoothie! :(
            length: length of gaussian window in unit of sigma!
        
            Returns
            ----------
        
            See also
            ----------
        
            Notes
            --------
            
    
            Last modified on Nov 2013
        '''
        ripples = []
        for item in self:
            if not hasattr(item, 'signal_filtered'):
                print 'Signal is not filtered yet!\n Filtering the signal with the default values...'
                item.filter(150, 250)
            if not hasattr(item, 'hilbertAbsolute'):
                print 'Hilbert transform is not calculated!\n Calculating Hilbert transform of the signal for' + item.tags['file'] + '...'
                item.hilbertTransform()
            hs_signal = numpy.zeros(self[0].hilbertAbsolute.size)  # averaging Hilbert signals!
        if sigma > 0:
            for item in self:
                item.hilbertAbsSmooth = numpy.convolve(item.hilbertAbsolute, scsig.gaussian(length * sigma, sigma), 'same')
        for i in range(len(self)):
            hs_signal = numpy.add(hs_signal, self[i].hilbertAbsSmooth)
        hs_signal /= float(len(self))       
        std = hs_signal.std()
        avg = hs_signal.mean()
        sig1 = numpy.where(hs_signal > (avg + std))[0]
        sig1_break = numpy.where(numpy.diff(sig1) > 1)[0]
        sig1_split = numpy.split(sig1, sig1_break + 1)
        t_chain = []
        t_start = 0.0
        t_end = 0.0
        peak = 0.0
        t_peak = 0.0
        for i in range(len(sig1_split)):
            if  hs_signal[sig1_split[i]].max() > avg + 3 * std:
                t_chain = numpy.where(hs_signal == hs_signal[sig1_split[i]].max())[0]
                t_start = self[0].timeAxis[sig1_split[i][0]]
                t_end = self[0].timeAxis[sig1_split[i][-1]]
                peak = hs_signal[sig1_split[i]].max()
                t_peak = self[0].timeAxis[numpy.intersect1d(sig1_split[i], t_chain)][0]
                ripples.append([t_start, t_end, peak, t_peak])
                
        ripples = numpy.array(ripples)
        self.ripples = ripples
        self.hs_signal = hs_signal
        print 'Ripple detection  is DONE!!!\n \
                Now you should see the attribute ripples with [t_start, t_end,peak value, t_peak] in each row for detected SWRs.'
        
    def badCSCs(self):
        '''
        A simple function to detect bad CSCs without visual inspection!
        Note: It only detects CSCs with unusual "quality". So if all CSCs are equally bad or good 
        it will not notice the difference!
        '''
        percent = numpy.array([])
        for i in range(len(self)):
            max_points = numpy.where(self[i].signal == self[i].signal.max())[0]
            percent = numpy.append(percent, max_points.size / float(self[i].signal.size) * 100)
        order_of_mag = numpy.log10(percent)
        thresh = 10 ** (order_of_mag.mean() + order_of_mag.std())
        for i in range(percent.size):
            if percent[i] > thresh:
                print self[i].tags['file']
        self.prcnt = percent

        
    
                    

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

