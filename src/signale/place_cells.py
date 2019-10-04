"""
signale.place_cells
===================

A module for place cells.
"""
__author__ = ("KT", "Moritz Dittmeyer")
__version__ = "4.3, September 2013"


# python modules

# other modules
import numpy
import matplotlib.pyplot as pl
import matplotlib.path as mpl_path
import NeuroTools.signals as NTsig

# custom made modules
import custom_plot

# package modules
from signals import signal
import tools



# colors
transred = '#FFA1A1'


###################################################### FUNCTIONS


def phasePlot(spikes, phases, places, traj, wHists=True, wTraj=True, fig=None, ax=None, labely=True, limx=None, color=[0, 0, 0]):

    if not fig:
        fig = pl.figure(figsize=(4.5, 3.5))
    if not ax:
        pos = [.25, .21, .65, .625]
        ax = fig.add_axes(pos)
    pos = ax.get_position().bounds

    color = numpy.array(color)

    ax.plot(places[:, 0], phases, '.', color=color)
    ax.plot(places[:, 0], phases + 360, '.', color=color)



    if limx is None or not limx:
        limx = [0, traj.trackLength]
    elif limx is True:
        dx = places[:, 0].max() - places[:, 0].min()
        dx *= .2
        limx = [max(0, places[:, 0].min() - dx), min(traj.trackLength, places[:, 0].max() + dx)]


    custom_plot.huebschMachen(ax)
    ax.set_xlabel('Position (' + traj.spaceUnit + ')')
    ax.set_xticks(numpy.linspace(numpy.round(limx[0], 1), numpy.round(limx[1], 1), 3))
    if labely:
        ax.set_ylabel('Spike phase (deg)')
        ax.set_yticks(range(0, 720 + 1, 180))
    else:
        ax.set_yticks(range(0, 720 + 1, 180))
        ax.set_yticklabels([])
    ax.set_xlim(limx[0], limx[1])
    ax.set_ylim(0, 720)


    if wHists:
        histCol = custom_plot.grau
        if color.any():
            histCol = color

        posAbove = list(pos)
        posAbove[1] = pos[1] + pos[3]
        posAbove[3] = .05
        posRight = list(pos)
        posRight[0] = pos[0] + pos[2]
        posRight[2] = .05
        axAbove = fig.add_axes(posAbove)
        axAbove.hist(places[:, 0], numpy.arange(0, traj.trackLength, .03), facecolor=histCol, edgecolor=custom_plot.grau * 2. / 3)
        custom_plot.allOff(axAbove)
        axAbove.set_xlim(limx[0], limx[1])

        axRight = fig.add_axes(posRight)

        axRight.hist(numpy.append(phases, phases + 360), numpy.arange(0, 720, 11.25), \
            orientation='horizontal', facecolor=histCol, edgecolor=custom_plot.grau * 2. / 3)
        axRight.set_ylim(0, 720)
        custom_plot.allOff(axRight)
        ax = axAbove

    if wTraj:
        posUppermost = list(posAbove)
        posUppermost[1] = posAbove[1] + posAbove[3]
        posUppermost[3] = .05
        axUppermost = fig.add_axes(posUppermost)
        axUppermost.plot(traj.places[:, 0], traj.places[:, 1], linewidth=1, color=numpy.ones(3) * .5)  # original/full trajectory

# #        if onlyRunning:
# #            axUppermost.plot(places[:, 0], places[:, 1], 'r.')
# #        else:
        axUppermost.plot(places[:, 0], places[:, 1], 'r.')

        axUppermost.set_xlim(limx[0], limx[1])
        axUppermost.set_ylim(traj.ylim + numpy.diff(traj.ylim) * .1 * [-1, 1])
        custom_plot.allOff(axUppermost)
        ax = axUppermost

    # plotTitle = stList.tags[zaehler]['file']+', '+cscName
    # ax.set_title(plotTitle, fontsize=fontsize-4)

    return fig, ax


###################################################### CLASSES


class placeCell_spikezug(signal, NTsig.SpikeTrain):
    """ Spike train class for a place modulated cell.

    It extends the NeuroTools SpikeTrain class.
    """

    def __init__(self, spike_times, t_start=None, t_stop=None, timeUnit='ms'):
        signal.__init__(self)
        NTsig.SpikeTrain.__init__(self, spike_times, t_start=None, t_stop=None)
        self.traj = None

    def getLeftAndRightwardSpikes(self, onlyRunning=False):

        inclPhases = True
        if not hasattr(self , 'spike_phases'):
            print 'WARNING: Calculating without spike phases! If spike phases are needed \
                spike phases have to be calculated first!'
            inclPhases = False

        if onlyRunning:
            st = numpy.copy(self.run_spikeTimes)
            rechts_spikeTimes = numpy.copy(self.run_spikeTimes)
            rechts_spikePlaces = numpy.copy(self.run_spikePlaces)
            links_spikeTimes = numpy.copy(self.run_spikeTimes)
            links_spikePlaces = numpy.copy(self.run_spikePlaces)

            if inclPhases:
                rechts_spikePhases = numpy.copy(self.run_spikePhases)
                links_spikePhases = numpy.copy(self.run_spikePhases)
        else:
            st = numpy.copy(self.spike_times)
            rechts_spikeTimes = numpy.copy(self.spike_times)
            rechts_spikePlaces = numpy.copy(self.spikePlaces)
            links_spikeTimes = numpy.copy(self.spike_times)
            links_spikePlaces = numpy.copy(self.spikePlaces)

            if inclPhases:
                rechts_spikePhases = numpy.copy(self.spike_phases)
                links_spikePhases = numpy.copy(self.spike_phases)

        for j, time in enumerate(st):
            if numpy.ma.count_masked(numpy.ma.masked_values(self.traj.rechts_times, \
                                     time, atol=self.traj.dt)) == 0:
                rechts_spikeTimes[j] = numpy.nan
                rechts_spikePlaces[j] = numpy.nan
                if inclPhases:
                    rechts_spikePhases[j] = numpy.nan
            if numpy.ma.count_masked(numpy.ma.masked_values(self.traj.links_times, \
                                     time, atol=self.traj.dt)) == 0:
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


    def getSpikePlaces(self, interp=True):
        """ Return the place at which a certain spike occured.
        """

        if interp:
            print "NOTE: Using linear interpolation to get the places."

        spikePlaces = []
        for spike in self.spike_times:
            if spike > self.traj.times[-1]:
                print 'WARNING: at time', spike, 'there was a spike after the\
                    end of the trajectory', self.traj.times[-1], '!'
                spikePlaces.append(numpy.ones(3) * 0)
            else:
                spikePlaces.extend(self.traj.getPlaceFromTime(spike, interp=interp))
        self.spikePlaces = numpy.array(spikePlaces)


    def getRunningSpikes(self, threshspeed=None):
        """ Get spikes that occured above certain running speed.
        """

        if not threshspeed:
            threshspeed = self.traj.threshspeed
        if threshspeed != self.traj.threshspeed:
            print 'NOTE: Reset threshspeed from', self.traj.threshspeed, 'to', \
                threshspeed, self.traj.spaceUnit + '/' + self.traj.timeUnit
            self.traj.threshspeed = threshspeed

            print 'NOTE: Now calculating running traj >', threshspeed, self.traj.spaceUnit + '/' + self.traj.timeUnit
            self.traj.getRunningTraj(threshspeed=threshspeed)


        if not hasattr(self.traj , 'run_times'):
            print 'NOTE: Now calculating running traj >', threshspeed, self.traj.spaceUnit + '/' + self.traj.timeUnit
            self.traj.getRunningTraj(threshspeed=threshspeed)


        # parameters
        inclPhases = True
        if not hasattr(self , 'spike_phases'):
            print 'WARNING: Calculating without spike phases!', \
                'If spike phases are needed caluculate them before!'
            inclPhases = False

        inclHeadDir = True
        if not hasattr(self , 'spike_headDirections'):
            print 'WARNING: Calculating without spike head directions!', \
                'If spike head directions are needed caluculate them before!'
            inclHeadDir = False

        # initialization
        self.run_spikeTimes = numpy.copy(self.spike_times)
        if inclPhases:
            self.run_spikePhases = numpy.copy(self.spike_phases)
        if inclHeadDir:
            self.run_spikeHeadDirections = numpy.copy(self.spike_headDirections)
        if not hasattr(self , 'spikePlaces'):
            print 'NOTE: Now calculating places of spikes.'
            self.getSpikePlaces()
        self.run_spikePlaces = numpy.copy(self.spikePlaces)


        for j, time in enumerate(self.spike_times):
            if numpy.ma.count_masked(numpy.ma.masked_values(self.traj.run_times, time, atol=self.traj.dt * 20)) == 0:
                self.run_spikeTimes[j] = numpy.nan
                self.run_spikePlaces[j] = numpy.nan
                if inclPhases:
                    self.run_spikePhases[j] = numpy.nan
                if inclHeadDir:
                    self.run_spikeHeadDirections[j] = numpy.nan

        # remove nans from arrays
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
            phases.append(lfp.hilbertPhase[tools.findNearest(lfp.timeAxis, spike)[0]])
        self.spike_phases = numpy.array(phases)

        return self.spike_phases


    def getSpikeHeadDirection(self, interp=True):
        """ Return the head direction at which a certain spike occured.
        """

        if interp:
            print "NOTE: Using interpolation to get head directions."

        spike_headDirections = []
        for spike in self.spike_times:
            if spike > self.traj.times[-1]:
                print 'WARNING: at time', spike, 'there was a spike after the trajectory was over!'
                spike_headDirections.append(0)
            else:
                spike_headDirections.extend(self.traj.getHeadDirectionFromTime(spike, interp=interp))
        self.spike_headDirections = numpy.array(spike_headDirections)


    def plotSpikesvsPlace(self, onlyRunning=False, showHeadDir=False, fig=None, ax=None):
        """ Plot spikes vs. trajectory.

        """

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        if not hasattr(self , 'spikePlaces'):
            self.getSpikePlaces()

        # plot trajectory
        ax.plot(self.traj.places[:, 0], self.traj.places[:, 1], linewidth=1, color=numpy.ones(3) * .5)  # original/full trajectory

        # plot spikes
        if not onlyRunning:
            ax.plot(self.spikePlaces[:, 0], self.spikePlaces[:, 1], color=transred, marker='.', linestyle='None')  # , label='all spikes')  # all spiketimes
        ax.plot(self.run_spikePlaces[:, 0], self.run_spikePlaces[:, 1], 'r.')  # , label='only running spikes') # only running-spikes

        # plot head direction
        if showHeadDir:
            hd = self.traj.getHeadDirection()
            ax.quiver(self.traj.places[:, 0], self.traj.places[:, 1], numpy.cos(hd), \
                numpy.sin(hd), pivot='mid', color='g', units='y')
            if not onlyRunning:
                ax.quiver(self.spikePlaces[:, 0], self.spikePlaces[:, 1], \
                    numpy.cos(self.spike_headDirections), \
                    numpy.sin(self.spike_headDirections), pivot='mid', color=[.1, .1, 1], units='dots')
            ax.quiver(self.run_spikePlaces[:, 0], self.run_spikePlaces[:, 1], \
                numpy.cos(self.run_spikeHeadDirections), \
                numpy.sin(self.run_spikeHeadDirections), pivot='mid', color=[0, 0, 1], units='dots')

        # huebsch machen
        ax.set_xlim(self.traj.xlim + numpy.diff(self.traj.xlim) * .05 * [-1, 1])
        ax.set_ylim(self.traj.ylim + numpy.diff(self.traj.ylim) * .1 * [-1, 1])

        ax.set_xlabel('x position (' + self.traj.spaceUnit + ')')
        ax.set_ylabel('y position (' + self.traj.spaceUnit + ')')

        return fig, ax


    def plotField(self, bin_size=.05, onlyRunning=False, fig=None, ax=None):
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
        track_x1 = numpy.arange(self.traj.xlim[0], self.traj.xlim[1], bin_size)
        track_x2 = numpy.arange(self.traj.xlim[0] + bin_size, self.traj.xlim[1] + bin_size, bin_size)
        track_y1 = numpy.arange(self.traj.ylim[0], self.traj.ylim[1], bin_size)
        track_y2 = numpy.arange(self.traj.ylim[0] + bin_size, self.traj.ylim[1] + bin_size, bin_size)

        # create polygones & count spikes
        spike_number = []
        for l, ySP in enumerate(track_y1):  #
            for j, xSP in enumerate(track_x1):
                pol = mpl_path.Path([[track_x1[j], track_y1[l]], [track_x2[j], track_y1[l]], \
                                    [track_x1[j], track_y2[l]], [track_x2[j], track_y2[l]]])  # Polygon erzeugen
                # pol.append([[track_x1[j],track_y1[l]], [track_x2[j],track_y1[l]], [track_x1[j],track_y2[l]], [track_x2[j],track_y2[l]]])
                spike_number.append(numpy.sum(pol.contains_points(spikes)))  # count number of spikes in polygon
        X, Y = numpy.meshgrid(track_x1, track_y1)
        spike_number = numpy.array(spike_number).reshape(len(track_y1), len(track_x1))  # reshape & array spike_number list
        # X,Y = numpy.meshgrid(numpy.arange(self.traj.xlim[0],self.traj.xlim[1]+bin_size, bin_size), numpy.arange(self.traj.ylim[0], self.traj.ylim[1]+bin_size, bin_size))


        ax.pcolor(X, Y, spike_number)
        
        ax.set_xlim(self.traj.xlim[0], self.traj.xlim[1] - bin_size)
        ax.set_ylim(self.traj.ylim[0], self.traj.ylim[1] - bin_size)

        ax.set_xlabel('x position (' + self.traj.spaceUnit + ')')
        ax.set_ylabel('y position (' + self.traj.spaceUnit + ')')

        return fig, ax



    def phasePlot(self, wHists=True, wTraj=True, fig=None, ax=None, labely=True, limx=False):

        fig, ax = phasePlot(self.spike_times, self.spike_phases, self.spikePlaces, self.traj, \
                            wHists=True, wTraj=True, fig=None, ax=None, labely=True, limx=False)

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
    
    
        
        
