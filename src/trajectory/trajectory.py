"""
trajectory.trajectory
=====================

Base module for spatial paths, i.e., trajectories.
"""

__author__ = ("KT", "Moritz Dittmeyer", "Benedikt Ludwig", "Sven Schoernich")
__version__ = "6.6, September 2013"

# python modules
import os, sys

# other modules
import numpy

import matplotlib.pyplot as pl
from matplotlib.collections import LineCollection
import matplotlib as mpl

# custom made modules
import custom_plot, signale


# package modules


###################################################### FUNCTIONS



def str2number(strs=[]):

    dummy = []
    if isinstance(strs, str):
        strs = [strs]

    for j, s in enumerate(strs):
        dummy.append(int(''.join(i for i in s if i.isdigit())))

    return dummy


def _change_folderStyle(dirName):
    """
    Change folder style from Windows to Linux and Python, respectively.
    """

    name = ''
    for d in dirName.split('\\'):
        name += d + '/'

    return name




def _getTimeUnitFactor(timeUnit, newUnit):
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


class trajectory(object):

    def __init__(self, traj, meta={}):

        # initialize meta data
        self.meta = meta
        self.dt = 0.1  # claim a dt (s)
        self.t_start = 0.0
        self.mazeType = ''
        self.timeUnit = 's'
        self.spaceUnit = 'm'
        self.view_eulerOffsets = [0, 0, 0]
        self.euler = -numpy.array(self.view_eulerOffsets)  # for storing the current euler angles
        self.euler %= 360
        self.threshspeed = 0  # minimum speed for running (m/s)


        # get real values if parsable from file
        if meta.has_key('dt'):
            self.dt = meta['dt']  # dt [s]
        if meta.has_key('t_start'):
            self.t_start = meta['t_start']  # t_start [s]
        if meta.has_key('time'):
            self.time = meta['time']  # time and date, when
                                                # recording was started,
                                                # just to provide some time stamp

        if meta.has_key('mazetype'):
            self.mazeType = meta['mazetype']

        self.times = []
        self.has_timeAxis = False
        if traj.shape[1] == 4:  # with time stamp?
            self.times = traj[:, 0]
            self.t_start = numpy.min(self.times)
            self.t_stop = numpy.max(self.times)
            self.places = traj[:, 1:]
            self.numPlaces = self.places.shape[0]
            self.has_timeAxis = True
        elif traj.shape[1] == 3:  # without time stamp?
            self.places = traj
            self.numPlaces = self.places.shape[0]
            self.t_stop = (self.numPlaces - 1) * self.dt + self.t_start
        else:
            sys.exit('Traj shape not implemented!')

        if self.dt > 0 and not self.times.__len__():  # calculate times if not provided
            self.times = numpy.arange(self.t_start, self.t_stop + self.dt, self.dt)  # gives the times, when the place is entered
            if self.times.shape[0] > self.places.shape[0]:
                self.times = self.times[:self.places.shape[0]]  # bad tweak, since sometimes t_stop
                                                                        # gets included by arange although it shouldn't!?

        self.turn(self.view_eulerOffsets[0])  # remove the yaw component of view_eulerOffsets


        # if dt was not provided, get it from times array
        if not meta.has_key('dt'):
            self.dt = numpy.mean(numpy.diff(self.times))
        # print "dt is", self.dt


    def __recalc_startstop(self):

        try: self.times
        except AttributeError:
            pass
        else:
            self.t_start = numpy.min(self.times)
            self.t_stop = numpy.max(self.times)


    def __changeUnits(self, newTimeUnit=None, newSpaceUnit=None):

        if newTimeUnit:
            factor = _getTimeUnitFactor(self.timeUnit, newUnit)

            # change the times
            self.times *= factor
            self.dt *= factor
            self.threshspeed /= factor
            self.timeUnit = newUnit

            self.__recalc_startstop()

        elif newSpaceUnit:
            print 'NOTE: Changing spaceUnit not implemented.'
        else:
            print 'NOTE: No units were changed.'



    def centerTraj(self):
        """
        Shifts the trajectory center to position (0,0).
        """
        self.getTrajDimensions()

        for i, p in enumerate(self.places):
            self.places[i, 0] += self.xWidth / 2
            self.places[i, 1] += self.yWidth / 2


    def cumPlaces(self, iStart=None, iStop=None):
        speed = self.getSpeed()[1]
        cumPlace = numpy.cumsum(speed[iStart:iStop]) * self.dt
        cumPlace = numpy.append([0], cumPlace)
        return cumPlace


    def getSpeed(self, thresh=None, vec=False, laps=False):
        """ Calculate running speed.

        thresh ... optional argument of minimum running speed

        Returns: time and speed numpy arrays
        """

        if laps:
            diffsTmp = []
            diffsTmp2 = []
            laps, lapTimes = self.getLaps()
            for item in lapTimes:
                lapStart = self.getIndexFromTime(item[0])
                lapEnd = self.getIndexFromTime(item[1])
                diffsTmp.append(numpy.diff(self.places[lapStart:lapEnd], axis=0)[:, 0:2])
                diffsTmp2.append(numpy.diff(self.times[lapStart:lapEnd]))

            diffs = []
            diffsTimes = []
            for index, item in enumerate(diffsTmp):
                for subitem in item:
                    diffs.append(subitem)
                diffsTimes.append(diffsTmp2[index])
            diffs = numpy.array(diffs)
            diffsTimes = numpy.array(diffsTimes)
        else:
            diffs = numpy.diff(self.places, axis=0)[:, 0:2]
            diffsTimes = numpy.diff(self.times)


        if vec:
            speed = diffs / diffsTimes
        else:
            speed = numpy.sqrt(diffs[:, 1] ** 2 + diffs[:, 0] ** 2) / diffsTimes  # [space units/s]

        speed_dummy = []
        reducedTimes = []
        if thresh and not vec:
            for i, s in enumerate(speed):
                if s >= thresh:
                    speed_dummy.append(s)
                    reducedTimes.append(self.times[i])
        else:
            speed_dummy = speed
            reducedTimes = self.times[:-1]
        return numpy.array(reducedTimes), numpy.array(speed_dummy)


    def getAcceleration(self, thresh=None, abs=False, vec=False, laps=False):
        """ Calculate acceleration.

        thresh ... optional argument to set a minimum

        Returns: time and acceleration numpy arrays
        """
        t, speed = self.getSpeed(vec=vec, laps=laps)

        diffs = numpy.diff(speed, axis=0)[:, 0:2]
        diffsTimes = numpy.diff(t)
        acceleration = diffs / diffsTimes  # [space units/s^2]

        if abs and not vec:
            acceleration = numpy.abs(acceleration)

        acceleration_dummy = []
        reducedTimes = []

        if thresh and not vec:
            for i, s in enumerate(acceleration):
                if numpy.abs(s) >= thresh:
                    acceleration_dummy.append(s)
                    reducedTimes.append(self.times[i])
        else:
            acceleration_dummy = acceleration
            reducedTimes = self.times[:-2]
        return numpy.array(reducedTimes), numpy.array(acceleration_dummy)


    def getJerk(self, thresh=None, abs=False, vec=False, laps=False):
        """ Calculate jerk, third derivative.

        thresh ... optional argument to set a minimum

        Returns: time and jerk numpy arrays
        """

        t, acceleration = self.getAcceleration(vec=vec, laps=laps)

        diffs = numpy.diff(acceleration, axis=0)[:, 0:2]
        diffsTimes = numpy.diff(acceleration)
        jerk = diffs / diffsTimes  # [space units/s^3]

        if abs and not vec:
            jerk = numpy.abs(jerk)

        jerk_dummy = []
        reducedTimes = []

        if thresh and not vec:
            for i, s in enumerate(jerk):
                if numpy.abs(s) >= thresh:
                    jerk_dummy.append(s)
                    reducedTimes.append(self.times[i])
        else:
            jerk_dummy = jerk
            reducedTimes = self.times[:-3]
        return numpy.array(reducedTimes), numpy.array(jerk_dummy)


    def getHeadDirection(self):
        """
        Calculates approximate head or gaze direction from trajectory.
        """

        diffs = numpy.diff(self.places, axis=0)
        lengths = numpy.sqrt(diffs[:, 0] ** 2 + diffs[:, 1] ** 2)
        angles = numpy.arccos(diffs[:, 0] / lengths) * numpy.sign(diffs[:, 0])

        indicesWithDivideByZero = numpy.nonzero(numpy.isnan(angles))[0]  # there will be NaNs

        # iterate through the angles and override NaN entries with the last real angle
        while indicesWithDivideByZero.shape[0] > 0:
            angles[indicesWithDivideByZero] = angles[indicesWithDivideByZero - 1]
            indicesWithDivideByZero = numpy.nonzero(numpy.isnan(angles))[0]

        self.headDirections = angles

        return self.headDirections


    def getIndexFromTime(self, time):
        """ Return index corresponding to the given time.

        Parameters:
        time ... a number or an 1D array

        Returns:
        index ... an int or array according to time
        """

        if not isinstance(time, list) and not isinstance(time, numpy.ndarray):
            # if necessary convert time value to list, i.e., array
            time = [time]

        index = []
        for t in time:
            index.append(signale.findNearest(self.times, t)[0])

        if index.__len__() == 1:  # convert index to int in case of one entry only
            index = index[0]
        else:
            index = numpy.array(index)

        return index


    def getPlaceFromTime(self, time, interp=False):
        """ Return place array corresponding to the given time array.
        """

        if not isinstance(time, list) and not isinstance(time, numpy.ndarray):
            # if necessary convert time value to list, i.e., array
            time = [time]

        place = []
        if interp:
            index1 = numpy.searchsorted(self.times, time, side='left') - 1  # get index of the time point that is just smaller than time
            index2 = index1 + 1
            t1 = self.times[index1]
            t2 = self.times[index2]
            tpercent = (time - t1) / (t2 - t1)

            place1 = self.places[index1]
            place2 = self.places[index2]
            place = tpercent * (place2 - place1) + place1
# #            print t1, t2, time, tpercent
# #            print place1, place2, place
# #            print
        else:
            index = self.getIndexFromTime(time)  # get index of the time point that is just smaller than time
            place.extend(self.places[[index]])

        return numpy.array(place)


    def getHeadDirectionFromTime(self, time, interp=False):
        """ Return head direction array corresponding to the given time array.
        """

        if not isinstance(time, list) and not isinstance(time, numpy.ndarray):
            # if necessary convert time value to list, i.e., array
            time = [time]

        headDirection = []
        if interp:
            index1 = numpy.searchsorted(self.times, time, side='left') - 1  # get index of the time point that is just smaller than time
            index2 = index1 + 1  # get index of the time point that is just bigger than time

            t1 = self.times[index1]
            t2 = self.times[index2]
            tpercent = (time - t1) / (t2 - t1)

            hd = self.getHeadDirection()
            hd1 = hd[index1]
            if index2 >= hd.size:
                index2 -= 1
            hd2 = hd[index2]
            headDirection = tpercent * numpy.absolute((numpy.exp(hd2 * 1j) - numpy.exp(hd1 * 1j))) + hd1
            # headDirection = tpercent * (hd2 - hd1) + hd1
        else:
            index = self.getIndexFromTime(time)  # get index of the time point that is just smaller than time
            headDirection.append(self.getHeadDirection()[index])

        return numpy.array(headDirection)


    def getRunningTraj(self, threshspeed=None, window_len=51):
        """ Get times & places, where the animal was running.

        threshspeed ... minimum running speed
        window_len ... window length for smooting

        """

        if not threshspeed:
            threshspeed = self.threshspeed
        if threshspeed != self.threshspeed:
            print 'NOTE: Reset threshspeed from', self.threshspeed, 'to', \
                threshspeed, self.spaceUnit + '/' + self.timeUnit
            self.threshspeed = threshspeed

        print "calculating trajectory with running speed >=", self.threshspeed

        inclHeadDir = True
        if not hasattr(self , 'headDirections'):
            print 'WARNING: Calculating without head directions! If head directions are needed \
                calculate them first!'
            inclHeadDir = False

        # get running speed and smooth it a bit
        speed_dummy = self.getSpeed()[1]  # signale.smooth(self.getSpeed()[1], window_len=window_len, window='hanning')

        indices = numpy.where(speed_dummy >= self.threshspeed)[0]
        indices += 1  # shift indices by one index
        if indices.size and indices[-1] >= speed_dummy.size:
            indices = indices[:-1]
        indices = numpy.append([0], indices)  # add entry at index zero to indices
        self.run_times = self.times[indices]
        self.run_places = self.places[indices]
        if inclHeadDir:
            self.run_headDirections = self.headDirections[indices]


    def getLaps(self):
        """
        Returns the number of laps in the trajectory.

        A lap is defined as the path/time between leaving a reward
        area and entering the next. NOTE, that the next reward area
        might be the same or a different reward area.
        """

        centers = self.rewardsPos
        radius = self.rewardsArea[0]
        places = self.places[:, 0:2]


        i = 0
        laps = 0
        inLap = False
        lapTimes = []

        while i < places.shape[0]:

            diffs = places[i] - centers
            distances = []
            for item in diffs:
                distances.append(numpy.sqrt(item[0] ** 2 + item[1] ** 2))
            minDistance = numpy.array(distances).min()

            if minDistance > radius and not inLap:
                laps += 1
                lapStart = self.times[i]
                inLap = True
            elif minDistance < radius and inLap:
                lapTimes.append([lapStart, self.times[i - 1]])
                inLap = False

            i += 1

        return laps, lapTimes


    def getTrajDimensions(self):
        """
        Detects coordinates with maximum/minimum x and y values
        and using this calculates the 'width' of the trajectory.
        """

        # self.orient()
        mini = self.places.argmin(0)
        maxi = self.places.argmax(0)

        xMin = self.places[mini[0], 0]
        xMax = self.places[maxi[0], 0]
        yMin = self.places[mini[1], 1]
        yMax = self.places[maxi[1], 1]

        self.xWidth = xMax - xMin
        self.yWidth = yMax - yMin

        self.xlim = numpy.array([xMin, xMax])
        self.ylim = numpy.array([yMin, yMax])

        return self.xWidth, self.yWidth


    def turn(self, yaw=0):
        """
        Turn trajectory by yaw degrees.
        """

        slope = numpy.tan(numpy.deg2rad(yaw))
        vecTrack = numpy.array([1 , slope])
        normalVecTrack = numpy.array([-vecTrack[1], vecTrack[0]])

        # do it for the places
        tmp = []
        for i, p in enumerate(self.places):
            tmp = p[0:2].copy()
            self.places[i, 0] = numpy.dot(tmp, vecTrack) / numpy.sqrt(numpy.dot(vecTrack, vecTrack))
            self.places[i, 1] = numpy.dot(tmp, normalVecTrack) / numpy.sqrt(numpy.dot(vecTrack, vecTrack))

        # do it again for the rewardsPos
        if hasattr(self, 'rewardsPos'):
            tmp = []
            for i, p in enumerate(self.rewardsPos):
                tmp = p[0:2].copy()
                self.rewardsPos[i, 0] = numpy.dot(tmp, vecTrack) / numpy.sqrt(numpy.dot(vecTrack, vecTrack))
                self.rewardsPos[i, 1] = numpy.dot(tmp, normalVecTrack) / numpy.sqrt(numpy.dot(vecTrack, vecTrack))

        self.getTrajDimensions()
        self.euler[0] += yaw


    def plot(self, fig=None, ax=None, offset=2, language='e', chic=False, marker=None):

        if not fig:
            fig = pl.figure(figsize=(8, 6))
        if not ax:
            ax = fig.add_axes([.2, .15, .81, .75])
        axcb = None

        if not marker:
            line_segments = LineCollection([[x, self.places[i + 1 + offset, 0:2]] \
                                for i, x in enumerate(self.places[offset:-(1 + offset), 0:2])], \
                                linestyles='solid', linewidths=mpl.rcParams['lines.linewidth'] / 2.)
            line_segments.set_array(self.times)
            ax.add_collection(line_segments)
        else:
            ax.plot(self.places[offset:, 0], self.places[offset:, 1], marker)

        ax.plot(self.places[0 + offset, 0], self.places[0 + offset, 1], 'o')  # start point
        ax.plot(self.places[-2, 0], self.places[-2, 1], 'd')  # end point


        # huebsch machen
        custom_plot.huebschMachen(ax)
        if not chic:
            ax.set_title(self.mazeType, fontsize=custom_plot.fontsize - 4)
        if language == 'd':  # in german
            ax.set_xlabel('x-Position (' + self.spaceUnit + ')')
            ax.set_ylabel('y-Position (' + self.spaceUnit + ')')
        else:  # in english
            ax.set_xlabel('x position (' + self.spaceUnit + ')')
            ax.set_ylabel('y position (' + self.spaceUnit + ')')
        dx = numpy.round((int(self.t_stop) / 60 * 60 - int(self.t_start) / 60 * 60) / 4.)
        if not dx:
            dx = 60
        xticks = numpy.arange(round(self.t_start), round(self.t_stop) + 1, dx)

        # colorbar
        if not marker:
            axcb = fig.colorbar(line_segments)
            axcb.set_ticks(xticks)
            if language == 'd':  # in german
                axcb.set_label('Zeit (' + self.timeUnit + ')')
            else:  # in english
                axcb.set_label('Time (' + self.timeUnit + ')')

        if not fig:
            pl.show()

        return fig, ax, axcb


    def plotCumPlaces(self, fig=None, ax=None, text=''):

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        cumPlaces = self.cumPlaces()
        ax.plot(self.times, cumPlaces)

        ax.set_xlabel('Time (' + self.timeUnit + ')')
        ax.set_ylabel('Path length (' + self.spaceUnit + ')')
        ax.text(self.times.max(), cumPlaces.max(), text)

        pl.show()

        return fig, ax


    def plotSpeed(self, thresh=None, fig=None, ax=None):

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)
        # ax = fig.add_subplot(211)

        reducedTimes, speed = self.getSpeed(thresh=thresh)
        ax.plot(reducedTimes, speed)

        ax.set_xlabel('Time (' + self.timeUnit + ')')
        ax.set_ylabel('Speed (' + self.spaceUnit + '/' + self.timeUnit + ')')
        ax.set_xlim(reducedTimes[0], reducedTimes[-1])

        # #--##

# #        ax = fig.add_subplot(212)
# #
# #        reducedTimes = self.getSpeed(thresh=thresh)[0]
# #        #reducedTimes = signale.smooth(self.getSpeed(thresh=thresh)[0], window_len=11, window='hanning')
# #        speed=signale.smooth(self.getSpeed(thresh=thresh)[1], window_len=11, window='hanning')
# #
# #        ax.plot(reducedTimes, speed)
# #
# #        ax.set_xlabel('Time ('+self.timeUnit+')')
# #        ax.set_ylabel('Smoothed-Speed ('+self.spaceUnit+'/'+self.timeUnit+')')
# #        ax.set_xlim(reducedTimes[0], reducedTimes[-1])
        pl.show()

        return fig, ax


    def plotSpeedvsPlace(self, thresh=None, fig=None, ax=None):

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        time, speed = self.getSpeed(thresh)
        places = self.getPlaceFromTime(time)

        offset = 2
        line_segments = LineCollection([[x, places[i + 1 + offset, 0:2]] for i, x in enumerate(places[offset:-2, 0:2])],
                               linestyles='solid')
        speed = numpy.minimum(speed, speed.mean() + 3.*speed.std())  # cut away outliers for plotting
        line_segments.set_array(speed)
        ax.add_collection(line_segments)

        fig = pl.gcf()
        axcb = fig.colorbar(line_segments)
        ax.plot(places[0 + offset, 0], places[0 + offset, 1], 'o')  # start point
        ax.plot(places[-2, 0], places[-2, 1], 'd')  # end point

        axcb.set_label('Speed (' + self.spaceUnit + '/' + self.timeUnit + ')')
        ax.set_title(self.mazeType)
        ax.set_xlabel('x position (' + self.spaceUnit + ')')
        ax.set_ylabel('y position (' + self.spaceUnit + ')')
        pl.show()

        return fig, ax


    def purge(self, numItems):
        """ Cut away numItems.
        """
        self.times = self.times[numItems:]
        self.places = self.places[numItems:]


    def removeValues(self, value=0.):
        """ Remove values due to artifacts (e.g., signal loss)
        """
        # for i in numpy.where(self.places[:,:2]==0)[0]:
            # self.places[i] = numpy.nan

        self.times = numpy.delete(self.times, numpy.where(self.places[:, :2] == value)[0], 0)
        self.places = numpy.delete(self.places, numpy.where(self.places[:, :2] == value)[0], 0)


    def recalcTimesAxis(self, places, times, time_before=8):
        """ Recalculate times.

        Find (average) smallest time offset to places in trajectory and
        recalculate the time axis of the trajectory accordingly.

        time_before default 8 (empirically tested)
        """
        if times.__len__():
            time_diffs = numpy.array([])
            for i, t in enumerate(times):
                index1 = self.getIndexFromTime(t) + 1
                index2 = self.getIndexFromTime(t - time_before)
                index_min = signale.findNearest_vec(self.places[index2:index1], places[i])[0] + index2
                time_diffs = numpy.append(time_diffs, t - self.times[index_min])
            self.time_offset(time_diffs.mean())
            # print "time differences:", time_diffs
            print "=> shifted time axis by" , time_diffs.mean(), ", std is" , time_diffs.std()
            self.__recalc_startstop()


    def time_slice(self, t0=None, t1=None):
        """ Return a trajectory object sliced by time.

        t0 ... start time for slice
        t1 ... end time for slice

        Returns: sliced trajectory
        """

        if not t0:
            t0 = self.times[0]

        if not t1:
            t1 = self.times[-1]

        # get indices corresponding to t0 and t1, respectively
        index0 = self.getIndexFromTime(t0)
        index1 = self.getIndexFromTime(t1)

        if index0 == index1:
            sys.exit('WARNING: Slice not possible!')

        sliced_traj_array = numpy.hstack([self.times.reshape(self.times.size, 1)[index0:index1, :], \
                                          self.places[index0:index1, :]])
        traj_type = type(self)
        sliced_traj = traj_type(sliced_traj_array, self.meta)

        return sliced_traj


    def time_offset(self, offset=0.):
        self.times += offset
        self.__recalc_startstop()




class linearMazeTrajectory(trajectory):

    def __init__(self, places, meta={}, initialized=False):
        if not initialized:
            trajectory.__init__(self, places, meta)  # initialized=True if trajectory.__init__() was already called
        self.oriented = False

        if not meta.has_key('expType'):
            print 'NOTE: No experiment type provided!'
            print '  => put some values as if experiment was done in real environment'
            meta['trackLength'] = 1.5
            meta['expType'] = 'real'
            meta['yaw'] = 0

        if meta.has_key('yaw'):
            self.yaw = meta['yaw'] % 360
            self.euler[0] += self.yaw
            self.euler %= 360

        self.orient()  # rotate trajectory to remove yaw
        if meta.has_key('trackLength') and meta.has_key('expType') and meta['expType'] == 'real':
            # self.orient()
            self.removeValues(0.)  # remove lines to zero
            self.getTrajDimensions()
            self.trackLength = self.xWidth
            self.trackWidth = self.yWidth
            for i, p in enumerate(self.places):
                self.places[i, 0:2] -= [self.xlim[0], self.ylim[0]]  # shift to get smallest place to 0,0
                self.places[i, 0:2] /= self.xWidth  # normalize it to trackLength
                self.places[i, 0:2] *= meta['trackLength']  # scale it to trackLength

        self.getTrajDimensions()
        self.trackLength = self.xWidth
        self.trackWidth = self.yWidth


    def getComponentsFromTime(self, time, interp=0):
        """
        Returns the x and y componets corresponding to the given time.
        """

        try: self.components
        except AttributeError:
            self.components = numpy.array([self.getXComponents(), self.getYComponents()]).T

        components = numpy.array([0, 0])
        if interp:
            index1 = numpy.searchsorted(self.times, time, side='right') - 1  # get index of the time point that is just smaller than time
            index2 = numpy.searchsorted(self.times, time, side='right')  # get index of the time point that is just bigger than time

            print 'NOTE: to do'
        else:
            index = numpy.searchsorted(self.times, time, side='right') - 1  # get index of the time point that is just smaller than time
            if abs(self.times[index] - time) > abs(self.times[index + 1] - time):
                index += 1
            components = self.components[index]

        return components


    def getLeftAndRightwardRuns(self, onlyRunning=True):

        if onlyRunning:
            if not hasattr(self , 'run_places'):
                print 'NOTE: Need to calculate running trajectory.'
                self.getRunningTraj()

            places = self.run_places
            times = self.run_times
        else:
            places = self.places
            times = self.times

        # x-Richtung bestimmen
        xdir = numpy.diff(signale.smooth(places[:, 0], window_len=11, window='hanning'), axis=0)

        self.rechts_places = numpy.copy(places)
        self.rechts_times = numpy.copy(times)
        self.links_places = numpy.copy(places)
        self.links_times = numpy.copy(times)

        indices = numpy.where(xdir <= 0)[0]
        self.rechts_places[indices] = numpy.nan
        self.rechts_times[indices] = numpy.nan

        indices = numpy.where(xdir > 0)[0]
        self.links_places[indices] = numpy.nan
        self.links_times[indices] = numpy.nan



    def getXComponents(self):
        xComp = numpy.array([])
        if hasattr(self, 'euler'):
            slope = -numpy.tan(numpy.deg2rad(self.euler[0]))
            vecTrack = numpy.array([1 , slope])
            for p in self.places:
                xComp = numpy.append(xComp, numpy.dot(p[0:2], vecTrack))
        else:
            vecTrack = self.rewardsPos[1] - self.rewardsPos[0]
            for p in self.places:
                xComp = numpy.append(xComp, numpy.dot(p[0:2] - self.rewardsPos[0], vecTrack))

        return xComp / numpy.sqrt(numpy.dot(vecTrack, vecTrack))


    def getYComponents(self):
        yComp = numpy.array([])
        if hasattr(self, 'euler'):
            slope = -numpy.tan(numpy.deg2rad(self.euler[0]))
            vecTrack = numpy.array([1 , slope])
            normalVecTrack = numpy.array([-vecTrack[1], vecTrack[0]])
            for p in self.places:
                yComp = numpy.append(yComp, numpy.dot(p[0:2], normalVecTrack))
        else:
            vecTrack = self.rewardsPos[1] - self.rewardsPos[0]
            normalVecTrack = numpy.array([-vecTrack[1], vecTrack[0]])
            for p in self.places:
                yComp = numpy.append(yComp, numpy.dot(p[0:2] - self.rewardsPos[0], normalVecTrack))

        return yComp / numpy.sqrt(numpy.dot(vecTrack, vecTrack))



    def orient(self):
        """
        Orients the trajectory to yaw = 0 deg by projecting it
        via a dot product with the slope vector.
        """
        if not self.oriented:
            self.turn(-self.euler[0])

            # rotate to positive axes, if necessary
            if self.yaw > 90 and self.yaw <= 270:
                self.places[:, 0:2] *= -1

            self.getTrajDimensions()
            self.trackWidth = self.yWidth
            self.oriented = True
        else:
            print "NOTE: Track already oriented."


    def plot(self, fig=None, ax=None, offset=2, language='e', chic=False):

        fig, ax, axcb = trajectory.plot(self, fig, ax, offset, language, chic)

        # adjust axes
        pos = [.15, .65, .65, .3]
        ax.set_position(pos)
        poscb = list(axcb.ax.get_position().bounds)
        poscb[0] = .825
        poscb[1] = pos[1]
        poscb[3] = pos[3]
        axcb.ax.set_position(poscb)


        # huebsch machen
        xoffset = self.xWidth * .1
        yoffset = self.yWidth

        xmin = -xoffset
        xmax = self.xWidth + xoffset
        ymin = self.ylim[0] - yoffset
        ymax = self.ylim[1] + yoffset

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)


        # show yaw
        if hasattr(self, 'yaw'):
            ax_inset = fig.add_axes([.73, .85, .07, .07], polar=True)
            x = numpy.deg2rad(self.yaw)
            patch = ax_inset.fill([x - numpy.pi, x - .85 * numpy.pi, x, x + .85 * numpy.pi, x + numpy.pi], [.75, 1, 1, 1, .75], facecolor='b', edgecolor='b')
            # custom_plot.drop_shadow_patches(ax_inset, patch[0])
            ax_inset.set_yticks([])
            ax_inset.set_xticks([])
            ax_inset.set_xticks(numpy.arange(0, 1, .125) * 2 * numpy.pi)
            ax_inset.set_xticklabels([])
            # ax_inset.set_xticklabels(numpy.int_(numpy.arange(0, 1, .25)*360), fontsize=custom_plot.fontsize/2)
            ax_inset.spines['polar'].set_linestyle('dotted')
            ax_inset.spines['polar'].set_linewidth(.5)

        # put additional axes
        pos1 = list(pos)
        pos1[1] = .11
        pos1[2] = .32
        pos1[3] = .35
        ax1 = fig.add_axes(pos1)
        pos2 = list(pos1)
        pos2[0] += pos2[2] + .16
        ax2 = fig.add_axes(pos2)

        xComponents = self.getXComponents()
        ax1.plot(self.times, xComponents, 'k-')

        if chic:
            yComponents = self.cumPlaces()
        else:
            yComponents = self.getYComponents()
        ax2.plot(self.times, yComponents, 'k-')

        # huebsch machen
        custom_plot.huebschMachen(ax1)
        custom_plot.huebschMachen(ax2)
        if language == 'd':  # in german
            ax1.set_xlabel('Zeit (' + self.timeUnit + ')')
            ax2.set_xlabel('Zeit (' + self.timeUnit + ')')
            ax1.set_ylabel('x-Position (' + self.spaceUnit + ')')
            ax2.set_ylabel('y-Position (' + self.spaceUnit + ')')
        else:  # in english
            ax1.set_xlabel('Time (' + self.timeUnit + ')')
            ax2.set_xlabel('Time (' + self.timeUnit + ')')
            ax1.set_ylabel('x position (' + self.spaceUnit + ')')
            if chic:
                ax2.set_ylabel('Path length (' + self.spaceUnit + ')')
            else:
                ax2.set_ylabel('y position (' + self.spaceUnit + ')')

        minY = xComponents.min()
        maxY = xComponents.max()
        dY = 1.
        ax1.set_yticks(numpy.arange(round(minY), round(maxY) + dY, dY))
        ax1.set_ylim(minY - dY / 10, maxY + dY / 10)
        dx = numpy.round((int(self.t_stop) / 60 * 60 - int(self.t_start) / 60 * 60) / 4.)
        if not dx:
            dx = 60
        xticks = numpy.arange(round(self.t_start), round(self.t_stop) + 1, dx)
        ax1.set_xticks(xticks)

        minY = yComponents.min()
        maxY = yComponents.max()
        dY = .1
        if chic:
            dY = int((maxY - minY) / 5)
        ax2.set_yticks(numpy.arange(round(minY, 1), round(maxY, 1) + dY, dY))
        ax2.set_ylim(minY - dY / 10, maxY + dY / 10)
        ax2.set_xticks(xticks)

        # customize plot further?
        if chic:
            custom_plot.allOff(ax)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')
            pos1[3] = .45
            pos2[3] = .45
            ax1.set_position(pos1)
            ax2.set_position(pos2)

            # add scale bars
            custom_plot.add_scalebar(ax, matchx=True, matchy=True, labelExtra=' m', \
                loc=3, borderpad=-1., sep=5)

        return fig, ax, axcb


    def plotComponents(self, plotInOne=0, plotAgainstEachOther=0, color='k', lineStyle='-'):
        """
        Plot x and y components of the linearMaze trajectory.

        plotInOne ... plot x and y components into one graph [default 0]
        plotAgainstEachOther ... plot x and y components against each other [default 0]
        color ... color of the plot [default 'k']
        lineStyle ... line style of the plot [default '-']

        Returns: nothing
        """

        ax = None
        ax2 = None

        if not plotInOne:
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.set_position([.15, .1, .8, .8])
            ax.plot(self.times, self.getXComponents(), color + lineStyle)
            ax.set_xlabel('Time (' + self.timeUnit + ')')
            ax.set_ylabel('x position (' + self.spaceUnit + ')')
            ax.set_xlim(self.times[0], self.times[-1])

            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.set_position([.15, .1, .8, .8])
            ax.plot(self.times, self.getYComponents(), color + lineStyle)
            ax.set_xlabel('Time (' + self.timeUnit + ')]')
            ax.set_ylabel('y position (' + self.spaceUnit + ')')
            ax.set_xlim(self.times[0], self.times[-1])
        else:
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.set_position([.15, .1, .7, .8])
            ax2 = ax.twinx()

            ax.plot(self.times, self.getXComponents(), 'k' + lineStyle)
            ax.set_xlabel('Time (' + self.timeUnit + ')')
            ax.set_ylabel('x position (' + self.spaceUnit + ')', color='k')
            ax.set_xlim(self.times[0], self.times[-1])
            for tl in ax.get_yticklabels():
                tl.set_color('k')

            ax2.plot(self.times, self.getYComponents(), 'k' + lineStyle, alpha=.5)
            ax2.set_ylabel('y position (' + self.spaceUnit + ')', color=[.5, .5, .5])
            for tl in ax2.get_yticklabels():
                tl.set_color([.5, .5, .5])

        if plotAgainstEachOther:
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.set_position([.15, .1, .8, .8])
            ax.plot(self.getXComponents(), self.getYComponents(), color + lineStyle)
            ax.set_xlabel('x position (' + self.spaceUnit + ')')
            ax.set_ylabel('y position (' + self.spaceUnit + ')')

        pl.show()

        return ax, ax2

    def plotSpeedvsComponents(self, thresh=None, avg=False, fig=None, ax=None, text='', color='k', lineStyle='.'):

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        time, speed = self.getSpeed(thresh)
        x = []
        y = []
        for t in time:
            dummy = self.getComponentsFromTime(t)
            x.append(dummy[0].tolist())
            y.append(dummy[1].tolist())
        x = numpy.array(x)
        y = numpy.array(y)

        if avg:
            bins = numpy.linspace(0, self.trackLength, 15)
            inds = numpy.digitize(x, bins)
            set_inds = numpy.unique(inds)
            average = numpy.zeros(set_inds.size)
            err = numpy.zeros(set_inds.size)
            for i, s in enumerate(set_inds):
                average[i] = speed[inds == s].mean()
                err[i] = speed[inds == s].std() / numpy.sqrt(numpy.sum(inds == s))
            custom_plot.avgPlot(ax, bins, average)  # , err=err
            ax.text(bins[average.argmax()], average.max(), text)
        else:
            ax.plot(x, speed, color + lineStyle)
        ax.set_xlabel('x position (' + self.spaceUnit + ')')
        ax.set_ylabel('Speed (' + self.spaceUnit + '/' + self.timeUnit + ')')
        ax.set_xlim(x.min(), x.max())

        if not avg:
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.plot(y, speed, color + lineStyle)
            ax.set_xlabel('y position (' + self.spaceUnit + ')')
            ax.set_ylabel('Speed (' + self.spaceUnit + '/' + self.timeUnit + ')')
            ax.set_xlim(y.min(), y.max())

        pl.show()

        return fig, ax


class paramsTrajectory(object):

    def __init__(self, times=[], parameters={}, meta={}):

        # initialize meta data
        self.meta = meta
        self.t_start = 0.0
        self.mazeType = ''
        self.timeUnit = 's'

        # get real values if parsable from file
        if meta.has_key('dt'):
            self.dt = meta['dt']  # dt [s]
        if meta.has_key('t_start'):
            self.t_start = meta['t_start']  # t_start [s]
        if meta.has_key('time'):
            self.time = meta['time']  # time and date, when
                                                # recording was started,
                                                # just to provide some time stamp
        if meta.has_key('mazetype'):
            self.mazeType = meta['mazetype']

        self.times = times
        self.parameters = parameters


    def getParameter(self, key):
        """ Returns an array for a parameter from the parameters dictionary."""

        param_array = []
        for p in self.parameters:
            param_array.append(p[key])

        return param_array
