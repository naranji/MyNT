'''
Created on Nov 27, 2013

@author: chenani
'''
# other modules
import os
import numpy as np
import matplotlib.pyplot as pl
import pickle as pkl
import signale
#import matplotlib.gridspec as gridspec
# custom made modules
def load(filename):
        '''
        This function loads ephys files into working memory.
        '''
        pkl.load(filename)
        print "%s has been loaded." % filename
# package modules



class ephys:
    '''
    A class that fuses spikeLists and CSC's
    '''


    def __init__(self, spikezugList, NeuralynxCSCList):
        '''
        
        '''
        self.spikes = spikezugList
        self.LFPs = NeuralynxCSCList
        self.timeAxis = NeuralynxCSCList[0].timeAxis
    def coincidence(self, csc_no=-1, window=300):
        '''
        This function detects ripple episode containing at least one spike!
        
        Parameters
        ----------
        csc_no : This is the csc object from LFPs, which its ripples gonna be used to detect coincidence.
        window : This is the size of gaussian window which will by ripple_recorder function if the selected csc has no ripple attribute yet.
        
        Returns
        ----------
        overlap : This the array containing [unit#,spike#,ripple#] on each row for every coinciding ripple/spike!
        ripples : Ripples detected/copied as an attribute to the ephys object.
        csc : link to the csc which is selected as a argument of function for coincidence detection purposes.
        See also
        ----------
        signale.cscs.NeuralynxCSC.ripple_recorder, numpy.where
        
        Notes
        ----------
        

        '''
        
        if csc_no == -1:
            ripples = self.LFPs.ripples
        elif (csc_no >= 0 and csc_no < len(self.LFPs)):
            if hasattr(self.LFPs, 'ripples'):
                ripples = self.LFPs[csc_no].ripples
            else:
                self.LFPs[csc_no].ripple_recorder(window)
                ripples = self.LFPs[csc_no].ripples
            ccdc = []
            for k in range(len(self.spikes)):
                for i in range(self.spikes[k].spike_times.size):
                    less = np.where(ripples[:, 0] < self.spikes[k].spike_times[i])[0]  # Ripples started before spike!    
                    more = np.where(ripples[:, 1] > self.spikes[k].spike_times[i])[0]  # Ripples ended after spike!
                    ripp_idx = np.intersect1d(less, more)  # index of the ripple covering the spike!
                    if ripp_idx.size:
                        ccdc.append([k, i, ripp_idx])
            ccdc = np.array(ccdc)
            ccdc = ccdc[np.argsort(ccdc[:, 2])]  # sorting w.r.t ripple index  
            self.overlap = ccdc
            self.ripples = ripples
            self.csc = self.LFPs[csc_no]
            print "The result is accessible at overlap attr with the format: [cell#,spike#,ripple#]"
        elif csc_no == 's':
            print 'Using surrogate data...'
            self.spikes.burstDetector()
            self.spikes.MuaDetector()
            mua = self.spikes.MultiUnitFiring
            ccdc = []
            for k in range(len(self.spikes)):
                for i in range(self.spikes[k].spike_times.size):
                    less = np.where(mua[:, 0] < self.spikes[k].spike_times[i])[0]  # Ripples started before spike!    
                    more = np.where(mua[:, 1] > self.spikes[k].spike_times[i])[0]  # Ripples ended after spike!
                    mua_idx = np.intersect1d(less, more)  # index of the ripple covering the spike!
                    if mua_idx.size:
                        ccdc.append([k, i, mua_idx])
            ccdc = np.array(ccdc)
            ccdc = ccdc[np.argsort(ccdc[:, 2])]  # sorting w.r.t ripple index  
            self.overlap = ccdc
        else:
            print 'cscIndex out of range!'     
           
        
        
    def MUA(self):
        '''
        This function detects ripple episodes containing spikes from at least two different cells.
        
        Parameters
        ----------
        There is no parameters, this function basically uses the "overlap" attr. if its already calculated by coincidence() finction.
        Otherwise it will call the coincidence() first and then does its job!!!
        
        Returns
        ----------
        rippleSpikeContent : 2D array with ripple number on the first column and number of spikes it contains on the second!
        MultiUnitActivity : Vector of arrays with each array containing spike information during a ripple. The format is [cell#,spike#,ripple#]
        
        See also
        ----------
        coincidence, numpy.where
        
        Notes
        ----------
         It takes the result of coincidence function as input.
        '''
        if not hasattr(self, 'overlap'):
            print 'There is no overlap detected between spikes and ripples, calculating using default arguments...'
            self.coincidence()
        mua = []  # Multi Unit Activity
        bursts = []  # This keeps the firing activity(>=2spikes) during ripples
        repetition = []  # Contains the number of spikes in each ripple!
        ripp_set = set(self.overlap[:, 2])  # Keeps the tags of ripples containing MUA! 
        
        for item in ripp_set:
            repetition.append([item, np.where(self.overlap[:, 2] == item)[0].size])
        repetition = np.array(repetition)
        self.rippleSpikeContent = repetition
        
        for item in repetition[:, 0]:
            activity_block = self.overlap[np.where(self.overlap[:, 2] == item)]
            if not activity_block[:, 0].mean() == activity_block[0, 0]:
                mua.append(activity_block)
        self.MultiUnitActivity = mua
        
    def sequencing(self):
        '''
        This function sequences the events detected by MUA() function
        
        Parameters
        ----------
        Returns
        ----------
        sequence_1st_spike : 
        sequence_median :
        
        See also
        ----------
        coincidence, MUA
        
        Notes
        ----------
         It takes the result of MUA function as input.
        '''
        fst_spk_sq = []
        median_sq = []  # Sequences based on first spike times!
        for item in self.MultiUnitActivity:
            cells = set(item[:, 0])
            sequence = []
            for jtem in cells:
                row_idx = np.where(item[:, 0] == jtem)[0]
                spk_times = self.spikes[jtem].spike_times[item[row_idx, 1]]  # This returns the spike times of a specific cell.
                sequence.append([jtem, spk_times.min(), np.median(spk_times), spk_times.size, item[0, 2]]) #No use for the last 2! Why did i put them there?
            sequence = np.array(sequence)
            sequence = sequence[np.argsort(sequence[:, 1])]
            fst_spk_sq.append(sequence)
            sequence = sequence[np.argsort(sequence[:, 2])]
            median_sq.append(sequence)
        self.sequence_1st_spike = fst_spk_sq
        self.sequence_median = median_sq
    def ripp_sq_plot(self, sequence_no=0, sequencing_method='median',fig = None, ax1 = None,ax2 = None, ax3 = None):
        '''
        This function plots the raster plot of a sequence during a specific ripple indicated by a pair [csc_no,ripp_no]
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        # ingredients :)
        if not fig:
            fig = pl.figure()
        if not ax1:
            ax1 = fig.add_subplot(311)
        if not ax2:
            ax2 = fig.add_subplot(312)
        if not ax3:
            ax3 = fig.add_subplot(313)
        spk_slice = []
        # Finding the desired ripple and assign it to link ripp.
        if sequencing_method == 'median':
            sqn = self.sequence_median[sequence_no]
            ripp_idx = sqn[0, -1]
            cell_sq = sqn[:, 0]
            ripp = self.ripples[ripp_idx] 
        elif sequencing_method == '1st':
            sqn = self.sequence_1st[sequence_no]
            ripp_idx = sqn[0, -1]
            cell_sq = sqn[:, 0]
            ripp = self.ripples[ripp_idx]
        # copy the time slice synced with ripple from the spike trains in the order of sequence! 
        for item in cell_sq:
            spk_slice.append(self.spikes[item].time_slice(ripp[0], ripp[1]))
        # Raster plots
        for k in range(len(spk_slice)):
            band = 1 / float(len(spk_slice))
            for i in range(len(spk_slice[k])):
                ax1.plot([spk_slice[k].spike_times[i], spk_slice[k].spike_times[i]], [k - 0.2, k + 0.2], 'k', linewidth=2)
        ax1.set_ylim(-1, len(spk_slice) - 0.5)
        ax1.set_xlim(ripp[0], ripp[1])
        ax1.set_yticks(range(len(spk_slice)))
        ax1.set_ylabel('Cells') 
        # Ripple plot
        t_i = np.where(self.timeAxis == ripp[0])[0][0]  # index of ripple onset
        t_f = np.where(self.timeAxis == ripp[1])[0][0]  # index of ripple end
        ax2.plot(self.timeAxis[t_i:t_f], self.csc.signal_filtered[t_i:t_f], 'r')
        ax2.set_ylim(self.csc.signal[t_i:t_f].min(), self.csc.signal[t_i:t_f].max())
        ax2.set_yticks([])
        # Raw signal
        ax3.plot(self.timeAxis[t_i:t_f]-self.timeAxis[t_i], self.csc.signal[t_i:t_f], 'k')
        ax3.set_yticks([])
        fig.subplots_adjust(hspace=0.0)
        return fig,ax1,ax2,ax3
    def rippleSpikeContent(self,method='fixed_window'):
        '''
        A function to find spike times during ripples!
        Parameters
        ----------
        
        
        Returns
        ----------
        ripple_MUA: List consist of vectors containing spike times for different spike trains during given ripple.
        
        See also
        ----------
        
        Notes
        ----------
        '''
        if not hasattr(self.LFPs, 'ripples'):
            print 'Ripples have not detected yet!!! \n Lets do that first...'
            self.LFPs.rippleDetect(sigma=10)
        else:
            if method=='fixed_window':
                rippPeak = np.array([item[3] for item in self.LFPs.ripples])
                ripple_MUA = []
                for ii in range(rippPeak.size):
                    rippSpikeTimes = []
                    for zug in self.spikes:
                        spk_idx = np.intersect1d(np.where(zug.spike_times < np.int64(rippPeak[ii]) + 50)[0],np.where(zug.spike_times > np.int64(rippPeak[ii]) - 50)[0])
                        if spk_idx.size > 0:
                            spkTms = zug.spike_times[spk_idx]
                        else:
                            spkTms = np.array([])
                        rippSpikeTimes.append(spkTms)
                    ripple_MUA.append(rippSpikeTimes)
                self.ripple_MUA = ripple_MUA
                
            if method=='dynamic_window':
                rippStart = np.array([item[0] for item in self.LFPs.ripples])
                rippEnd =   np.array([item[1] for item in self.LFPs.ripples])
                ripple_MUA = []
                for ii in range(rippEnd.size):
                    rippSpikeTimes = []
                    for zug in self.spikes:
                        spk_idx = np.intersect1d(np.where(zug.spike_times < rippEnd[ii])[0],np.where(zug.spike_times > rippStart[ii])[0])
                        if spk_idx.size > 0:
                            spkTms = zug.spike_times[spk_idx]
                        else:
                            spkTms = np.array([])
                        rippSpikeTimes.append(spkTms)
                    ripple_MUA.append(rippSpikeTimes)
                self.ripple_MUA = ripple_MUA
            
            
    def sequence_raster(self, sequence_no=0, sequencing_method='median',fig = None, ax1 = None,placeFields=np.array([]),placeSort = False,ClrMap=[]):
        '''
        This function plots the raster plot of a sequence!
        Parameters
        ----------
        
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        # ingredients :)
        if not fig:
            fig = pl.figure()
        if not ax1:
            ax1 = fig.add_subplot(111)
        spk_slice = []
        times = []
        
        for item in self.MultiUnitActivity[sequence_no]:
            cell = item[0]
            time_idx = item[1]
            times.append(self.spikes[cell].spike_times[time_idx])
        times = np.array(times)
        
        
        if sequencing_method == 'median':
            sqn = self.sequence_median[sequence_no]
            cell_sq = np.int0(sqn[:, 0])
            if placeFields.size > 0 and placeSort:
                cell_sq = np.int0([item for item in placeFields if item  in  cell_sq])
            if placeFields.size > 0 and not placeSort:
                cell_sq = np.int0([item for item in cell_sq if item  in  placeFields])
        elif sequencing_method == '1st':
            sqn = self.sequence_1st_spike[sequence_no]
            cell_sq = np.int0(sqn[:, 0])
            if placeFields.size > 0 and placeSort:
                cell_sq = np.int0([item for item in placeFields if item  in  cell_sq])
            if placeFields.size > 0 and not placeSort:
                cell_sq = np.int0([item for item in cell_sq if item  in  placeFields])
        # copy the time slice synced with ripple from the spike trains in the order of sequence! 
        for item in cell_sq:
            spk_slice.append(self.spikes[item].time_slice(times.min(), times.max()))
        # Raster plots
        if len(ClrMap):
            for k,cell in enumerate(spk_slice):
                band = 1 / float(len(spk_slice))
                for i in range(len(spk_slice[k])):
                    ax1.plot([spk_slice[k].spike_times[i], spk_slice[k].spike_times[i]], [k - 0.2, k + 0.2], color=ClrMap[cell_sq[k]], linewidth=2)
        else:
            for k,cell in enumerate(spk_slice):
                band = 1 / float(len(spk_slice))
                for i in range(len(spk_slice[k])):
                    ax1.plot([spk_slice[k].spike_times[i], spk_slice[k].spike_times[i]], [k - 0.2, k + 0.2], 'k', linewidth=2)

        ax1.set_ylim(-1, len(spk_slice) - 0.5)
        ax1.set_xlim(times.min()-1, times.max()+1)
        ax1.set_yticks(range(len(spk_slice)))
        ax1.set_ylabel('Cells') 
        if placeSort:
            ax1.set_title('Place Sorted')
        return fig,ax1
    def save(self, filename):
        '''
        This function saves the ephys object to a .ephys file!
        Parameters
        ----------
        filename:
        
        Returns
        ----------
        
        See also
        ----------
        
        Notes
        ----------
        '''
        pkl.dump(self, open(filename + ".ephys", 'wb'), pkl.HIGHEST_PROTOCOL)
        print "data has been saved to %s.ephys" % filename
    
    
