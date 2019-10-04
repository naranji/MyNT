"""
A module for spatial paths, i.e., trajectories.
"""

__author__ = ("KT", "Benedikt Ludwig", "Sven Schoernich")
__version__ = "3.31, June 2012"

# python modules
import sys

# other modules
import numpy

import matplotlib.pyplot as pl
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
import matplotlib as mpl

# custom made modules
import custom_plot


###################################################### FUNCTIONS


def decisionLine(n, h=None, p0=.5, p1=.75, alpha=.05, beta=.01):
    """
    h ... hypothesis type, 0 for H0 and 1 for H1

    Returns:
    decision line
    """
    if not h in [0, 1]:
        print 'Calculate for H0 or H1?'
        return 0

    x = numpy.log((1 - p0) / (1 - p1))
    nenner = numpy.log(p1 / p0) + x

    b = x / nenner
    if h:
        a = numpy.log((1 - beta) / alpha) / nenner
    else:
        a = numpy.log(beta / (1 - alpha)) / nenner

    return b * n + a


def find_nearest(a, v):
    """ Find the entry nearest to a value in an 1D array.

    Paramters:
    a ... 1D array
    v ... value to find nearest entry for

    Returns:
    a tuple (index, a[index])
    """
    index = (numpy.abs(a - v)).argmin()
    return (index, a[index])


def find_nearest_vec(a, v):
    """ Find the entry nearest to a vector value in an 1D array.

    Paramters:
    a ... 1D array
    v ... value to find nearest entry for

    Returns:
    a tuple (index, a[index])
    """
    diff = a - v
    norms = numpy.array([])
    for d in diff:
        norms = numpy.append(norms, numpy.linalg.norm(d))
    index = norms.argmin()
    return (index, a[index])


def load_trajectory(fileName, showHeader=False):

    traj = numpy.loadtxt(fileName)
    if traj.shape.__len__() == 1:
        traj = traj.reshape(1, traj.shape[-1])
    meta = _read_metadata(fileName, showHeader=showHeader)

    if meta['mazetype'] == 'linearMaze':
        return linearMazeTrajectory(traj, meta)
    elif meta['mazetype'].find('yDecision') + 1:
        return decisionMazeTrajectory(traj, meta)
    else:
        return trajectory(traj, meta)


def load_collisionTrajectory(fileName):

    traj = numpy.loadtxt(fileName)
    if traj.shape.__len__() == 1:
        traj = traj.reshape(1, traj.shape[-1])

    meta = _read_metadata(fileName)
    if meta['mazetype'] == 'linearMaze':
        return linearMazeCollisionTrajectory(traj, meta)
    else:
        return collisionTrajectory(traj, meta)


def load_rewardTrajectory(fileName):

    traj = numpy.loadtxt(fileName)
    if traj.shape.__len__() == 1:
        traj = traj.reshape(1, traj.shape[-1])
    meta = _read_metadata(fileName)

    if meta['mazetype'] == 'linearMaze':
        return linearMazeRewardTrajectory(traj, meta)
    elif meta['mazetype'].find('yDecision') + 1:
        return decisionMazeRewardTrajectory(traj, meta)
    else:
        return rewardTrajectory(traj, meta)


def load_rewardsPos(fileName):
    """
    For loading reward related information for a maze coming from Blender.
    """

    # initialize
    rewardsPosNames = []
    rewardsPosNamesAdditionals = {}
    rewardsPos = []
    initialPos = []
    rewardsArea = []
    rewardsAreaType = ''

    meta = _read_metadata(fileName)
    for key in meta:
        if key == 'rewardsArea':
            exec 'dummy=' + str(meta[key])
            if not isinstance(dummy, list):
                dummy = [dummy]
            rewardsArea.extend(dummy)
        elif key == 'rewardsRadius':
            rewardsArea = [meta[key]]
            print 'NOTE: loading reward positions with deprecated header style using rewardsRadius'
        elif key == 'rewardsAreaType':
            rewardsAreaType = meta['rewardsAreaType']
        elif key.endswith('add'):
            exec 'dummy=' + meta[key]
            rewardsPosNamesAdditionals[key.split('add')[0]] = dummy
        elif key == 'initialPos':
            exec 'dummy=' + meta[key]
            initialPos.extend(dummy)
        else:
            rewardsPosNames.append(key)
            exec 'dummy=' + meta[key]
            rewardsPos.append(dummy)
    if not len(initialPos):
        initialPos = rewardsPos[0]
        print 'NOTE: loaded obsolete maze style without initialPos.'

    return rewardsPosNames, rewardsPosNamesAdditionals, rewardsPos, rewardsArea, rewardsAreaType, initialPos


def load_userCFG(fileName):
    """ Load user config file
    """
    cfg = _read_metadata(fileName)
    return cfg


def save_trajectory(fileName, traj, meta={}):
    file = open(fileName, 'w')

    # write header by scanning through the meta dictionary and saving items appropriately
    _save_metadata(file, meta)

    if meta.has_key('trackEndsPos'):
        print 'NOTE: writing data files with deprecated header style using trackEndsPos'
    if meta.has_key('trackEndsRadius'):
        print 'NOTE: writing data files with deprecated header style using trackEndsRadius'
# #    if meta.has_key('start_time'):
# #        file.write( '# start_time=\''+meta['start_time']+'\'\n')
# #    if meta.has_key('end_time'):
# #        file.write( '# end_time=\''+meta['end_time']+'\'\n')
# #    if meta.has_key('dt'):
# #        file.write( '# dt=%f \n' % meta['dt'])                # dt [s]
# #    if meta.has_key('t_start'):
# #        file.write( '# t_start=%f \n' % meta['t_start'])
# #    if meta.has_key('mazetype'):
# #        file.write( '# mazetype=\''+meta['mazetype']+'\'\n')
# #        if meta['mazetype']=='linearMaze':
# #            if meta.has_key('yaw'):
# #                file.write( '# yaw=%f \n' % meta['yaw'])
# #            if meta.has_key('trackEndsPos'):
# #                print 'NOTE: writing data files with old header style using trackEndsPos'
# #                file.write( '# trackEndsPos=\''+str(meta['trackEndsPos'])+'\'\n')
# #            if meta.has_key('trackEndsRadius'):
# #                print 'NOTE: writing data files with old header style using trackEndsRadius'
# #                file.write( '# trackEndsRadius=%f \n' % meta['trackEndsRadius'])
# #    if meta.has_key('rewardsPos'):
# #        file.write( '# rewardsPos=\''+str(meta['rewardsPos'])+'\'\n')
# #    if meta.has_key('rewardsRadius'):
# #        file.write( '# rewardsRadius=%f \n' % meta['rewardsRadius'])


    # write actual data
    numpy.savetxt(file, traj, fmt='%s')
    file.close()

def save_userCFG(fileName, userCFG):
    """ Load user config file
    """
    file = open(fileName, 'w')
    _save_metadata(file, userCFG)
    file.close()

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


def _read_metadata(fileName, showHeader=False):
    """
    Read the informations that may be contained in the header of
    the trajectory object, if saved in a text file
    """
    metadata = {}
    cmd = ''
    f = open(fileName, 'r')
    for line in f.readlines():
        if line[0] == '#':
            cmd += line[1:].strip() + ';'
            if showHeader:
                print line.strip()
        else:
            break
    f.close()
    exec cmd in None, metadata

    if showHeader:
        print ''

    return metadata


def _save_metadata(file, meta={}):
    """ Save metadata.

    NOTE: Expects an file opened for reading by the calling function.
    """
    for key in meta:
        if isinstance(meta[key], float):
            file.write('# ' + str(key) + '=%f \n' % meta[key])  # str(key) is just for savety, in case key is not a string
        elif isinstance(meta[key], int):
            file.write('# ' + str(key) + '=%i \n' % meta[key])
        elif isinstance(meta[key], str) or isinstance(meta[key], unicode):
            file.write('# ' + str(key) + '=\'' + meta[key] + '\'\n')
        elif isinstance(meta[key], list):
            file.write('# ' + str(key) + '=\'' + str(meta[key]) + '\'\n')


###################################################### CLASSES


class trajectory(object):

    def __init__(self, traj, meta):

        # initialize meta data
        self.meta = meta
        self.dt = 0.0
        self.t_start = 0.0
        self.mazeType = ''
        self.timeUnit = 's'
        self.spaceUnit = 'm'

        # get real values if parseable from file
        if meta.has_key('dt'):
            self.dt = meta['dt']  # dt [s]
        if meta.has_key('t_start'):
            self.t_start = meta['t_start']  # t_start [s]
        if meta.has_key('time'):
            self.time = meta['time']  # time and date, when recording was started, just to provide some time stamp
        if meta.has_key('mazetype'):
            self.mazeType = meta['mazetype']
        
        self.view_eulerOffsets = [0, 0, 0]
        if meta.has_key('view_eulerOffsets'):
            exec 'self.view_eulerOffsets=' + meta['view_eulerOffsets']
        self.euler = -numpy.array(self.view_eulerOffsets)  # for storing the current euler angles
        self.euler %= 360
        
        if meta.has_key('rewardsPos'):  # the reward positions
            exec 'self.rewardsPos=' + meta['rewardsPos']
            self.rewardsPos = numpy.array(self.rewardsPos)
        elif meta.has_key('trackEndsPos'):  # actually not the track ends but the reward positions
            print 'NOTE: data files contain old header style with trackEndsPos'
            exec 'self.rewardsPos=' + meta['trackEndsPos']
            self.rewardsPos = numpy.array(self.rewardsPos)

        if meta.has_key('rewardsRadius'):  # the radius of the reward positions
            self.rewardsArea = numpy.array([meta['rewardsRadius']])
            print 'NOTE: data file contains old header style with rewardsRadius'
        elif meta.has_key('rewardsArea'):  # the extent of the reward positions
            exec 'self.rewardsArea=' + meta['rewardsArea']
            if not isinstance(self.rewardsArea, list):
                self.rewardsArea = [self.rewardsArea]
            self.rewardsArea = numpy.array(self.rewardsArea)
        elif meta.has_key('trackEndsRadius'):  # actually the radius of the reward positions
            print 'NOTE: data file contains old header style with trackEndsRadius'
            self.rewardsArea = numpy.array([meta['trackEndsRadius']])

        if meta.has_key('rewardsAreaType'):  # the extent of the reward positions
            self.rewardsAreaType = meta['rewardsAreaType']
        elif len(self.rewardsArea) == 1:
            self.rewardsAreaType = 'circular'
            print 'NOTE: no rewardsAreaType given, assuming circular!'
        elif len(self.rewardsArea) == 2:
            self.rewardsAreaType = 'rectangular'
            print 'NOTE: no rewardsAreaType given, assuming rectangular!'

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
    
    def __recalc_startstop(self):

        try: self.times
        except AttributeError:
            pass
        else:
            self.t_start = numpy.min(self.times)
            self.t_stop = numpy.max(self.times)


    def centerTrack(self):
        """
        Shifts the track (trajectory) center to position (0,0).
        """
        self.getTrajDimensions()

        for i, p in enumerate(self.places):
            self.places[i, 0] += self.xWidth / 2
            self.places[i, 1] += self.yWidth / 2

        for i, p in enumerate(self.rewardsPos):
            self.rewardsPos[i, 0] += self.xWidth / 2
            self.rewardsPos[i, 1] += self.yWidth / 2


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

        return self.xWidth, self.yWidth


    def cumPlaces(self, iStart=None, iStop=None):
        speed = self.getSpeed()[1]
        cumPlace = numpy.cumsum(speed[iStart:iStop]) * self.dt
        return cumPlace


    def getSpeed(self, thresh=None):
        """ Calculate running speed.

        thresh ... optional argument of minimum running speed

        Returns: time and speed numpy arrays
        """

        diffs = numpy.diff(self.places, axis=0)
        speed = numpy.sqrt(diffs[:, 1] ** 2 + diffs[:, 0] ** 2) / self.dt  # [space units/s]
        speed_dummy = []
        reducedTimes = []
        if thresh:
            for i, s in enumerate(speed):
                if s >= thresh:
                    speed_dummy.append(s)
                    reducedTimes.append(self.times[i])
        else:
            speed_dummy = speed
            reducedTimes = self.times[:-1]
        return numpy.array(reducedTimes), numpy.array(speed_dummy)

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

        return angles


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
            index.append(find_nearest(self.times, t)[0])
# #        index=numpy.searchsorted(self.times,time,side='right')-1    # get index of the time point that is just smaller than time
# #        if index+1<self.times.shape[0] and abs(self.times[index]-time) > abs(self.times[index+1]-time):
# #            index+=1

        if index.__len__() == 1:  # convert index to int in case of one entry only
            index = index[0]
        else:
            index = numpy.array(index)

        return index


    def getPlaceFromTime(self, time, interp=0):
        """ Return place array corresponding to the given time array.
        """

        if not isinstance(time, list) and not isinstance(time, numpy.ndarray):
            # if necessary convert time value to list, i.e., array
            time = [time]

        place = []
        if interp:
            index1 = numpy.searchsorted(self.times, time, side='right') - 1  # get index of the time point that is just smaller than time
            index2 = numpy.searchsorted(self.times, time, side='right')  # get index of the time point that is just bigger than time

            print 'NOTE: to do'
        else:
            # for t in time:
                # index=numpy.searchsorted(self.times, t, side='right')-1    # get index of the time point that is just smaller than time
                # if abs(self.times[index]-t) > abs(self.times[index+1]-t):
                #    index+=1
            #    place.append(self.places[index])
            index = self.getIndexFromTime(time)
            place.extend(self.places[[index]])
            place = numpy.array(place)
        
        return place


    def turn(self, yaw=0):
        """
        Turn track (trajectory) yaw degrees.
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
        tmp = []
        for i, p in enumerate(self.rewardsPos):
            tmp = p[0:2].copy()
            self.rewardsPos[i, 0] = numpy.dot(tmp, vecTrack) / numpy.sqrt(numpy.dot(vecTrack, vecTrack))
            self.rewardsPos[i, 1] = numpy.dot(tmp, normalVecTrack) / numpy.sqrt(numpy.dot(vecTrack, vecTrack))
    
        self.getTrajDimensions()
        self.euler[0] += yaw


    def plot(self, fig=None, ax=None, offset=2, language='e'):

        if not fig:
            fig = pl.figure(figsize=(8, 7))
        if not ax:
            ax = fig.add_subplot(111)

        line_segments = LineCollection([[x, self.places[i + 1 + offset, 0:2]] \
                            for i, x in enumerate(self.places[offset:-offset, 0:2])], \
                            linestyles='solid')
        line_segments.set_array(self.times)
        ax.add_collection(line_segments)

        axcb = fig.colorbar(line_segments)
        ax.plot(self.places[0 + offset, 0], self.places[0 + offset, 1], 'o')  # start point
        ax.plot(self.places[-2, 0], self.places[-2, 1], 'd')  # end point


        # huebsch machen
        custom_plot.huebschMachen(ax)
        ax.set_title(self.mazeType)
        if language == 'd':  # in german
            axcb.set_label('Zeit (' + self.timeUnit + ')')
            ax.set_xlabel('x-Position (' + self.spaceUnit + ')')
            ax.set_ylabel('y-Position (' + self.spaceUnit + ')')
        else:  # in english
            axcb.set_label('Time (' + self.timeUnit + ')')
            ax.set_xlabel('x position (' + self.spaceUnit + ')')
            ax.set_ylabel('y position (' + self.spaceUnit + ')')

        if not fig:
            pl.show()

        return fig, ax, axcb


    def plotCumPlaces(self):
        fig = pl.figure()
        ax = fig.add_subplot(111)

        cumPlace = self.cumPlaces()
        ax.plot(self.times[0:-1], cumPlace)

        ax.set_xlabel('Time (' + self.timeUnit + ')')
        ax.set_ylabel('Cummulative path length (' + self.spaceUnit + ')')
        pl.show()

    def plotSpeed(self, thresh=None):
        fig = pl.figure()
        ax = fig.add_subplot(111)

        reducedTimes, speed = self.getSpeed(thresh=thresh)
        ax.plot(reducedTimes, speed)

        ax.set_xlabel('Time (' + self.timeUnit + ')')
        ax.set_ylabel('Speed (' + self.spaceUnit + '/' + self.timeUnit + ')')
        ax.set_xlim(reducedTimes[0], reducedTimes[-1])
        pl.show()

    def plotSpeedvsPlace(self, thresh=None):
        time, speed = self.getSpeed(thresh)
        places = self.getPlaceFromTime(time)

        fig = pl.figure()
        ax = fig.add_subplot(111)

        offset = 2
        line_segments = LineCollection([[x, places[i + 1 + offset, 0:2]] for i, x in enumerate(places[offset:-2, 0:2])],
                               linestyles='solid')
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

    def recalcTimesAxis(self, places, times, time_before=7):
        """ Recalculate times.

        Find (average) smallest time offset to places in trajectory and
        recalculate the time axis of the trajectory accordingly.
        """
        if times.__len__():
            time_diffs = numpy.array([])
            for i, t in enumerate(times):
                index1 = self.getIndexFromTime(t) + 1
                index2 = self.getIndexFromTime(t - time_before)
                index_min = find_nearest_vec(self.places[index2:index1], places[i])[0] + index2
                time_diffs = numpy.append(time_diffs, t - self.times[index_min])
            self.time_offset(time_diffs.mean())
            print "time differences:", time_diffs
            print "=> shifted time axis by" , time_diffs.mean()
            self.__recalc_startstop()

    def purge(self, numItems):
        """
        Cut away numItems.
        """
        self.times = self.times[numItems:]
        self.places = self.places[numItems:]

    def time_offset(self, offset=0.):
        self.times += offset
        self.__recalc_startstop()


class linearMazeTrajectory(trajectory):

    def __init__(self, places, meta, initialized=False):
        if not initialized:
            trajectory.__init__(self, places, meta)  # initialized=True if trajectory.__init__() was already called
        if meta.has_key('yaw'):
            self.yaw = meta['yaw'] % 360
            self.euler[0] += self.yaw
            self.euler %= 360  

        # self.orient()      # rotate trajectory to remove yaw

        self.getTrajDimensions()

        # get track length
        # it actually isn't the track length but the distance between the reward positions
        # expected to be at the track ends
        self.trackLength = numpy.diff(self.rewardsPos, axis=0)
        self.trackLength = numpy.sqrt(self.trackLength[0, 0] ** 2 + self.trackLength[0, 1] ** 2)
        self.trackLength -= 2 * self.rewardsArea[0]
        # NOTE: previous command assumes circular reward zone with radius rewardsArea[0]

        self.trackWidth = self.yWidth


    def getXComponents(self):
        xComp = numpy.array([])
        if hasattr(self, 'euler'):
# #            yaw = self.yaw - self.view_eulerOffsets[0]
# #            if self.oriented:
# #                yaw = 0
# #            slope = -numpy.tan(numpy.deg2rad(yaw))
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
# #            yaw = self.yaw - self.view_eulerOffsets[0]
# #            if self.oriented:
# #                yaw = 0
# #            slope = -numpy.tan(numpy.deg2rad(yaw))
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
        Orients the track (trajectory) to yaw = 0 deg by projecting it
        via a dot product with the slope vector.
        """
        self.turn(-self.euler[0])
        
        # rotate to positive axes, if necessary
        if self.yaw > 90 and self.yaw <= 270:
            self.places[:, 0:2] *= -1
            self.rewardsPos[:, 0:2] *= -1

        self.getTrajDimensions()
        # get track length
        # it actually isn't the track length but the distance between the reward positions
        # expected to be at the track ends
        self.trackLength = numpy.diff(self.rewardsPos, axis=0)
        self.trackLength = numpy.sqrt(self.trackLength[0, 0] ** 2 + self.trackLength[0, 1] ** 2)
        self.trackLength -= 2 * self.rewardsArea[0]
        # NOTE: previous command assumes circular reward zone with radius rewardsArea[0]

        self.trackWidth = self.yWidth

# #    def orient(self):
# #        """
# #        Orients the track (trajectory) to yaw = 0 deg by projecting it
# #        via a dot product with the slope vector.
# #        """
# #        dummy0 = self.getXComponents()
# #        dummy1 = self.getYComponents()
# #        self.places[:, 0] = dummy0
# #        self.places[:, 1] = dummy1
# #        
# #
# #        # do it again for the rewardsPos
# #        yaw = self.yaw - self.view_eulerOffsets[0]
# #        if self.oriented:
# #                yaw = 0
# #        slope = -numpy.tan(numpy.deg2rad(yaw))
# #        vecTrack = numpy.array([1 , slope])
# #        normalVecTrack=numpy.array([-vecTrack[1], vecTrack[0]])
# #        tmp = []
# #        for i, p in enumerate(self.rewardsPos):
# #            tmp = p[0:2].copy()
# #            self.rewardsPos[i, 0] = numpy.dot(tmp, vecTrack)/numpy.sqrt(numpy.dot(vecTrack,vecTrack))
# #            self.rewardsPos[i, 1] = numpy.dot(tmp, normalVecTrack)/numpy.sqrt(numpy.dot(vecTrack,vecTrack))
# #
# #        # rotate to positive axes, if necessary
# #        if self.yaw > 90 and self.yaw < 270:
# #            self.places[:, 0:2] *= -1
# #            self.rewardsPos[:, 0:2] *= -1
# #
# #        self.getTrajDimensions()
# #        # get track length
# #        # it actually isn't the track length but the distance between the reward positions
# #        # expected to be at the track ends
# #        self.trackLength = numpy.diff(self.rewardsPos, axis=0)
# #        self.trackLength = numpy.sqrt(self.trackLength[0,0]**2+self.trackLength[0,1]**2)
# #        self.trackLength -= 2*self.rewardsArea[0]
# #        # NOTE: previous command assumes circular reward zone with radius rewardsArea[0]
# #
# #        self.trackWidth = self.yWidth
# #
# #        self.oriented = True


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

    def plot(self, fig=None, ax=None, offset=2, language='e'):

        fig, ax, axcb = trajectory.plot(self, fig, ax, offset, language)

        # adjust axes
        pos = [.15, .53, .65, .4]
        ax.set_position(pos)
        poscb = list(axcb.ax.get_position().bounds)
        poscb[0] = .825
        poscb[1] = pos[1]
        poscb[3] = pos[3]
        axcb.ax.set_position(poscb)


        # draw rewardsPos
        for rPos in self.rewardsPos:
            if self.rewardsAreaType == 'circular':
                circ = pl.Circle(rPos, radius=self.rewardsArea[0], color=custom_plot.grau2)
                ax.add_artist(circ)
            elif self.rewardsAreaType == 'rectangular':
                pass

        # huebsch machen
        xoffset = .5
        yoffset = .5

        xmin = -xoffset
        xmax = self.trackLength + 2 * self.rewardsArea[0] + xoffset
        ymin = -yoffset
        ymax = self.trackWidth + yoffset

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)


        # show yaw
        if hasattr(self, 'yaw'):
            # Convert to radians and subtract half the width
            # of a bar to center it.
            width = 0.3  # width of the bars (in radians)
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
        pos[1] = .1
        pos[2] = .32
        pos[3] = .3
        ax1 = fig.add_axes(pos)
        pos[0] += pos[2] + .16
        ax2 = fig.add_axes(pos)

        xComponents = self.getXComponents()
        yComponents = self.getYComponents()

        ax1.plot(self.times, xComponents, 'k-')
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
            ax2.set_ylabel('y position (' + self.spaceUnit + ')')

        minY = xComponents.min()
        maxY = xComponents.max()
        dY = 1.
        ax1.set_yticks(numpy.arange(round(minY), round(maxY) + dY, dY))
        ax1.set_ylim(minY - dY / 10, maxY + dY / 10)
        # xticks = numpy.arange(self.t_start, self.t_stop+1, numpy.round((self.t_stop-self.t_start)/5.))
        # ax1.set_xticks(xticks)

        minY = yComponents.min()
        maxY = yComponents.max()
        dY = .1
        ax2.set_yticks(numpy.arange(round(minY, 1), round(maxY, 1) + dY, dY))
        ax2.set_ylim(minY - dY / 10, maxY + dY / 10)
        # ax2.set_xticks(xticks)


    def plotComponents(self, plotInOne=0, plotAgainstEachOther=0, color='k', lineStyle='-'):
        """
        Plot x and y components of the linearMaze trajectory.

        plotInOne ... plot x and y components into one graph [default 0]
        plotAgainstEachOther ... plot x and y components against each other [default 0]
        color ... color of the plot [default 'k']
        lineStyle ... line style of the plot [default '-']

        Returns: nothing
        """

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

    def plotSpeedvsComponents(self, thresh=None, color='k', lineStyle='.'):
        time, speed = self.getSpeed()  # to do: including thresh
        x = self.getXComponents()[1:]
        y = self.getYComponents()[1:]

        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, speed, color + lineStyle)
        ax.set_xlabel('x position (' + self.spaceUnit + ')')
        ax.set_ylabel('Speed (' + self.spaceUnit + '/' + self.timeUnit + ')')
        ax.set_xlim(x.min(), x.max())

        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.plot(y, speed, color + lineStyle)
        ax.set_xlabel('y position (' + self.spaceUnit + ')')
        ax.set_ylabel('Speed (' + self.spaceUnit + '/' + self.timeUnit + ')')
        ax.set_xlim(y.min(), y.max())

        pl.show()

    def getLaps(self):
        """
        Returns the number of laps in the trajectory.

        A lap is defined as the path/time between leaving a reward
        area (track end) and entering the next. The next reward area
        might be the same or a different reward area.
        """
        x = self.getXComponents()
        i = 0
        laps = 0
        inLap = False
        lapTimes = []
        while i < x.shape[0]:
            if numpy.abs(x[i]) > 0 and numpy.abs(x[i]) < self.trackLength and not inLap:
                laps += 1
                lapStart = self.times[i]
                inLap = True
            elif inLap and (numpy.abs(x[i]) < 0 or numpy.abs(x[i]) > self.trackLength):
                lapTimes.append([lapStart, self.times[i - 1]])
                inLap = False
            i += 1
        return laps, lapTimes


class decisionMazeTrajectory(trajectory):

    def __init__(self, places, meta):
        trajectory.__init__(self, places, meta)
        if meta.has_key('yaw'):
            self.yaw = meta['yaw'] % 360
            self.euler[0] += self.yaw
            self.euler %= 360        
        
        self.wait_times = []
        if meta.has_key('wait_time'):
            self.wait_times.append(meta['wait_time'])
        if meta.has_key('wait_time_wrong'):
            self.wait_times.append(meta['wait_time_wrong'])
        self.wait_times = numpy.array(self.wait_times)


    def plot(self, times=[], fig=None, ax=None, offset=2, language='e'):

        if not fig:
            fig = pl.figure(figsize=(8, 7))
        if not ax:
            ax = fig.add_subplot(111)

        if not len(times):
            fig, ax, axcb = trajectory.plot(self, fig, ax, offset, language)
            return fig, ax, axcb


        indices1 = numpy.array([offset])
        indices2 = numpy.array([])
        indices1 = numpy.hstack((indices1, self.getIndexFromTime(times + self.wait_times[0]) + 1))
        indices2 = numpy.hstack((indices2, self.getIndexFromTime(times), self.times.size - 1))
        for i, idx1 in enumerate(indices1[:-1]):
            if indices1[i + 1] <= idx1:
                indices1[i + 1] = idx1 + 1
            if indices2[i] <= idx1:
                indices2[i] = idx1 + 1                
        print indices1, indices2
        print indices1 - indices2
        print times
        
        lines = []
        for i, index1 in enumerate(indices1):
            line = []
            for k in numpy.arange(index1, indices2[i]):
                line.append(tuple(self.places[k, 0:2].tolist()))
                line.append(tuple(self.places[k + 1, 0:2].tolist()))
            lines.append(tuple(line))

        # print len(lines), numpy.array(lines).shape, numpy.array(lines)[0]

        line_segments = LineCollection(lines, linestyles='solid')
        werte = numpy.arange(times.size + 2)
        line_segments.set_array(werte)
        ax.add_collection(line_segments)

        axcb = fig.colorbar(line_segments, boundaries=werte - .5)
        bins = werte.size
        if bins > 6:
            bins = 6
        werte = numpy.int_(numpy.linspace(werte[0], werte[-1] - 1, bins)) + 1
        axcb.set_ticks(werte - 1)
        axcb.set_ticklabels(werte)

        # huebsch machen
        custom_plot.huebschMachen(ax)
        ax.set_title(self.mazeType)
        if language == 'd':  # in german
            axcb.set_label('Lauf #')
            ax.set_xlabel('x-Position (' + self.spaceUnit + ')')
            ax.set_ylabel('y-Position (' + self.spaceUnit + ')')
        else:  # in english
            axcb.set_label('Run #')
            ax.set_xlabel('x position (' + self.spaceUnit + ')')
            ax.set_ylabel('y position (' + self.spaceUnit + ')')

        if not fig:
            pl.show()

        return fig, ax, axcb

# #    def purgeTrajectory(self, times):
# #        """ Remove parts during breaks after an rewardsPos has been entered,
# #            i.e., a decision has been made. The time is defined by wait_time. """
# #
# #        for t in times:
# #            index1 = find_nearest(self.times, t)[0]
# #            index2 = find_nearest(self.times, t+self.meta['wait_time'])[0]
# #            #self.times = numpy.delete(self.times, range(index1, index2+1))
# #            #self.places = numpy.delete(self.places, range(index1, index2+1), 0)
# #            self.times[index1:index2+1]=numpy.nan
# #            self.places[index1:index2+1]=numpy.nan
# #
# #
# #        self.purged=True


class collisionTrajectory(trajectory):

    def __init__(self, traj, meta):
        trajectory.__init__(self, traj, meta)
        # self.times=traj[:,0]

    def getICI(self, ICImin=-1):  # get Intercollision Interval
        """ Inter collision interval.
        """
        if ICImin == -1:
            ICImin = self.dt

        ICI = self.times[1:] - self.times[0:-1]
        self.reducedICI = numpy.array([])
        self.reducedPlaces = numpy.array([])
        self.reducedTimes = numpy.array([])
        for i, ic in enumerate(ICI):
            if ic >= ICImin:
                self.reducedICI = numpy.append(self.reducedICI, ic)
                self.reducedPlaces = numpy.append(self.reducedPlaces, self.places[i + 1])
                self.reducedTimes = numpy.append(self.reducedTimes, self.times[i + 1])
        return self.reducedICI

    def plotICI(self, ICImin=-1):
        fig = pl.figure()
        ax = fig.add_subplot(111)

        ICI = self.getICI(ICImin)
        ax.plot(ICI)

        ax.set_xlabel('Collision #')
        ax.set_ylabel('Intercollision interval [s]')
        pl.show()

    def getICD(self, traj, ICDmin=0, withTime=False):  # get Intercollision distance
        """Inter collision distance.
        """
        indices = numpy.array([])
        distances = numpy.array([])
        distances_time = numpy.array([])
        cumDistances = numpy.array([])
        for time in self.times:
            index = traj.getIndexFromTime(time)
            indices = numpy.append(indices, index)
        for i, index in enumerate(indices[0:-1]):
            if index < indices[i + 1]:  # just a sanity check
                d = traj.cumPlaces(index, indices[i + 1])
                if d[-1] >= ICDmin:
                    cumDistances = numpy.append(cumDistances, d)
                    distances = numpy.append(distances, d[-1])
                    if withTime:
                        distances_time = numpy.append(distances_time, self.times[i + 1])
        if withTime:
            return distances_time, distances, cumDistances
        else:
            return distances, cumDistances

    def plotICD(self, traj, ICDmin=0):
        [ICD_times, ICD, cumICD] = self.getICD(traj, ICDmin, withTime=True)

        fig = pl.figure()
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        # ax.plot(ICD_times, ICD, 'o')
        ax.plot(ICD)
        ax.set_xlabel('Collision #')
        # ax.set_xlim(self.t_start, self.t_stop)
        ax.set_xlim(0, ICD.size)
        
        ax2.plot(numpy.linspace(self.t_start, self.t_stop, cumICD.size), cumICD)

        ax2.set_xlabel('Time (' + self.timeUnit + ')')
        ax.set_ylabel('Intercollision distance (' + self.spaceUnit + ')')
        ax2.set_ylabel('Cumulated distance btw. collisions (' + self.spaceUnit + ')')
        ax2.set_xlim(self.t_start, self.t_stop)
        custom_plot.huebschMachen(ax)
        custom_plot.huebschMachen(ax2)
        pl.show()

    def purgeCollisions(self, dt=None):
        """Delete double entries (i.e., too close entries from collision trajectory)."""
        if not dt:
            dt = self.dt

        indices = numpy.where(numpy.diff(self.times) < dt)
        if indices[0].size:
            self.times = numpy.delete(self.times, indices)
            self.places = numpy.delete(self.places, indices, 0)

        self.numPlaces = self.places.shape[0]


class linearMazeCollisionTrajectory(collisionTrajectory, linearMazeTrajectory):

    def __init__(self, traj, meta):
        collisionTrajectory.__init__(self, traj, meta)
        linearMazeTrajectory.__init__(self, self.places, meta, initialized=True)
        # self.times=traj[:,0]


class rewardTrajectory(trajectory):

    def __init__(self, traj, meta):
        dummy = numpy.hstack((traj[:, 0:1], traj[:, 2:]))
        trajectory.__init__(self, dummy, meta)
        self.IDs = traj[:, 1]

    def getIRI(self):  # get Interreward interval
        "Inter reward interval"
        return self.times[1:] - self.times[0:-1]

    def plotIRI(self):
        fig = pl.figure()
        ax = fig.add_subplot(111)

        IRI = self.getIRI()
        ax.plot(IRI)

        ax.set_xlabel('Reward #')
        ax.set_ylabel('Inter reward interval [s]')
        pl.show()

    def getIRD(self, traj):  # get Interreward distance
        """Get the inter reward distance

        Parameters:
        traj ... the corresponding normal trajectory

        Returns:
        distances, cumDistances

        """
        indices = numpy.array([])
        distances = numpy.array([])
        cumDistances = numpy.array([])
        for time in self.times:
            index = traj.getIndexFromTime(time)
            indices = numpy.append(indices, index)
        for i, index in enumerate(indices[0:-1]):
            d = traj.cumPlaces(index, indices[i + 1])
            cumDistances = numpy.append(cumDistances, d)
            if d.size:
                distances = numpy.append(distances, d[-1])
        return distances, cumDistances

    def plotIRD(self, traj):
        """Plot the inter reward distance

        Parameters:
        traj ... the corresponding normal trajectory

        """
        [IRD, cumIRD] = self.getIRD(traj)

        fig = pl.figure()
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        ax.plot(IRD)
        ax2.plot(cumIRD)

        # draw track length
        if self.mazeType == 'linearMaze':
            ax.plot([0, IRD.shape[0]], [self.trackLength, self.trackLength])
            ax2.plot([0, cumIRD.shape[0]], [self.trackLength, self.trackLength])

        ax2.set_xlabel('Reward #')
        ax.set_ylabel('Inter reward distance (' + self.spaceUnit + ')')
        ax2.set_ylabel('Cummulated distance btw. rewards (' + self.spaceUnit + ')')
        pl.show()

    def purgeRewards(self, dt=None):
        """Delete double entries (i.e., too close entries from reward trajectory)."""
        if not dt:
            dt = self.dt

        indices = numpy.where(numpy.diff(self.times) < dt)
        if indices[0].size:
            self.times = numpy.delete(self.times, indices)
            self.places = numpy.delete(self.places, indices, 0)

        self.numPlaces = self.places.shape[0]


class linearMazeRewardTrajectory(rewardTrajectory, linearMazeTrajectory):

    def __init__(self, traj, meta):
        rewardTrajectory.__init__(self, traj, meta)
        linearMazeTrajectory.__init__(self, self.places, meta, initialized=True)
        # self.times=traj[:, 0]


class decisionMazeRewardTrajectory(rewardTrajectory):

    def __init__(self, traj, meta):
        # 1st entry of traj is the time
        # 2nd is the hotspot (stimulus pos, rewards pos) ID
        # from 4 we have the positions
        dummy = numpy.hstack((traj[:, 0:2], traj[:, 3:]))
        rewardTrajectory.__init__(self, dummy, meta)

        self.task = ''
        if meta.has_key('task'):
            self.task = meta['task']
        if meta.has_key('taskStimulus'):
            self.taskStimulus = meta['taskStimulus']
        if meta.has_key('stimuli'):
            exec 'self.stimuli=' + meta['stimuli']

        # get state of the stimuli
        # the third column is the state (binary) of the right arm
        self.leftStimulusIndex = traj[:, 2]
        self.rightStimulusIndex = abs(traj[:, 2] - 1)

    def getRightDecisions(self, display=False):

        rightDecisions = []
        for i, id in enumerate(self.IDs):
            if id:
                chosenColor = self.leftStimulusIndex[i]
            else:
                chosenColor = self.rightStimulusIndex[i]
            rightDecisions.append(int(self.taskStimulus == chosenColor))

        self.rightDecisions = numpy.array(rightDecisions)

        if display:
            fig = pl.figure(facecolor='w')
            ax = fig.add_axes([.15, .1, .7, .7])
            axSmall = fig.add_axes([.15, .85, .7, .07])

            bar_width = 0.8

            hist, bins = numpy.histogram(self.rightDecisions, [0, 1, 1.1])
            hist = numpy.float_(hist)
            bins = bins[:-1]
            bins -= bar_width / 2
            patches = ax.bar(bins, hist / hist.sum() * 100, width=bar_width, color=custom_plot.grau, edgecolor=custom_plot.grau)
            for p in patches:
                custom_plot.drop_shadow_patches(ax, p)
            ax.text(0, 5, str(int(hist[0])), horizontalalignment='center')
            ax.text(1, 5, str(int(hist[1])), horizontalalignment='center')
            ax.text(0, -7, 'wrong', horizontalalignment='center')
            ax.text(1, -7, 'right', horizontalalignment='center')

            custom_plot.huebschMachen(ax)
            ax.set_xticks([])
            ax.set_ylabel('Relative frequency (%)')
            ax.set_xlim(-.5, 1.5)
            ax.set_ylim(0, 100)
            ax.spines['bottom'].set_visible(False)

            patch = axSmall.axvspan(bins[0], bins[0] + bar_width, color=self.stimuli[abs(self.taskStimulus - 1)])
            custom_plot.drop_shadow_patches(ax, patch)
            patch = axSmall.axvspan(bins[1], bins[1] + bar_width, color=self.stimuli[self.taskStimulus])
            custom_plot.drop_shadow_patches(ax, patch)
            if self.task == 'visual':
                axSmall.text(bins[0] + bar_width / 2, 1.1, str(self.stimuli[abs(self.taskStimulus - 1)]), \
                        fontsize=mpl.rcParams['font.size'] - 4, horizontalalignment='center')
                axSmall.text(bins[1] + bar_width / 2, 1.1, str(self.stimuli[self.taskStimulus]), \
                        fontsize=mpl.rcParams['font.size'] - 4, horizontalalignment='center')
            custom_plot.allOff(axSmall)
            axSmall.set_xlim(-.5, 1.5)
            axSmall.text(1.5, .25, 'n=' + str(self.rightDecisions.size))

        return self.rightDecisions


    def sequentialAnalysis(self, display=False):

        if not hasattr(self, 'rightDecisions'):
            self.getRightDecisions()

        cumsum = self.rightDecisions.cumsum()
        n = numpy.arange(self.rightDecisions.size)
        alphas = [.05, .01, .001]
        h1_line = []
        h0_line = []
        h1True_index = []
        h0True_index = []
        h1True = []
        h0True = []
        for i, alpha in enumerate(alphas):
            h1_line.append(decisionLine(n, h=1, alpha=alpha))
            h0_line.append(decisionLine(n, h=0, alpha=alpha))

            # test if h1 gets true
            h1True_indices = numpy.where(h1_line[-1] - cumsum <= 0)[0]
            if h1True_indices.size:
                h1True_index.append(h1True_indices[0])
                h1True.append(n[h1True_indices[0]])
            else:
                h1True_index.append(-1)
                h1True.append(-1)

            # test if h0 gets true
            h0True_indices = numpy.where(h0_line[-1] - cumsum >= 0)[0]
            if h0True_indices.size:
                h0True_index.append(h0True_indices[0])
                h0True.append(n[h0True_indices[0]])
            else:
                h0True_index.append(-1)
                h0True.append(-1)

        if display:
            fig = pl.figure(facecolor='w')
            ax = fig.add_axes([.15, .12, .8, .68])
            axSmall = fig.add_axes([.15, .85, .8, .07])

            for i, alpha in enumerate(alphas):
                ax.plot(n, h1_line[i], color=custom_plot.colors[i], linewidth=3, alpha=.75, label=str(alpha))
                ax.plot(n, h0_line[i], color=custom_plot.colors[i], linewidth=3, alpha=.75)
            ax.plot(cumsum, 'bo-')
            for i, alpha in enumerate(alphas):
                if h1True_index[i] >= 0:
                    ax.plot(n[h1True_index[i]], cumsum[h1True_index[i]], 'ro')
                if h0True_index[i] >= 0:
                    ax.plot(n[h0True_index[i]], cumsum[h0True_index[i]], 'ro')
            custom_plot.huebschMachen(ax)
            ax.set_xlabel('Decision #')
            ax.set_ylabel('Right decisions')
            ax.set_xlim(0, self.rightDecisions.size)
            ax.legend(loc='upper left', frameon=False)

            if self.task == 'visual':
                axSmall.axvspan(0.25, .75, color=self.stimuli[abs(self.taskStimulus - 1)])
                axSmall.axvspan(1.25, 1.75, color=self.stimuli[self.taskStimulus])
                axSmall.text(0.5, 1.1, str(self.stimuli[abs(self.taskStimulus - 1)]), horizontalalignment='center')
                axSmall.text(1.5, 1.1, str(self.stimuli[self.taskStimulus]), horizontalalignment='center')
            axSmall.text(0.5, -.6, 'wrong', horizontalalignment='center')
            axSmall.text(1.5, -.6, 'right', horizontalalignment='center')
            custom_plot.allOff(axSmall)
            axSmall.set_xlim(0, 2)

        return numpy.array(h1True), numpy.array(h0True)



###################################################### traj_recorder


class traj_recorder(object):
    """ Base class for recording trajectories.
    
    Records arrays meant to be
    positions coordinates [x, y, z] and
    euler angles [yaw, pitch, roll].
    """

    def __init__(self, pos=[0., 0., 0.], eul=[0., 0., 0.], meta={}, realtimeplotting=0):
        # initialize meta data
        self.meta = meta
        self.dt = 0.0
        self.t_start = 0.0
        self.mazeType = ''

        # get real values if parseable
        if meta.has_key('dt'):
            self.dt = meta['dt']  # dt [s]
        if meta.has_key('t_start'):
            self.t_start = meta['t_start']  # t_start [s]
        if meta.has_key('start_time'):
            self.start_time = meta['start_time']  # time and date, when recording was started, just to provide some time stamp
        if meta.has_key('end_time'):
            self.end_time = meta['end_time']  # time and date, when recording was started, just to provide some time stamp
        if meta.has_key('mazetype'):
            self.mazeType = meta['mazetype']
        if self.mazeType == 'linearMaze':
            if meta.has_key('yaw'):
                self.yaw = meta['yaw']
            if meta.has_key('trackEndsPos'):
                exec 'self.trackEndsPos=' + meta['trackEndsPos']
                self.trackEndsPos = numpy.array(self.trackEndsPos)
            if meta.has_key('trackEndsRadius'):
                self.trackEndsRadius = meta['trackEndsRadius']

        # list of raw data
        self.posTraj = numpy.array(pos)
        self.eulTraj = numpy.array(eul)
        self.t = numpy.array(self.t_start)

        self.realtimeplotting = realtimeplotting
        self.headdirection = False

        # prepare figure for realtime plotting
        if self.realtimeplotting:
            pl.ion()
            self.fig = pl.figure()
            self.ax = self.fig.add_subplot(111)
            self.canvas = self.ax.figure.canvas
            self.ax.grid()

            # get line of data to be able to extend it later
            self.line = Line2D([self.posTraj[0]], [self.posTraj[1]])
            self.ax.add_line(self.line)

            # show initial head direction
            angle = numpy.pi / 2 - eul[0] * numpy.pi / 180
            pl.arrow(pos[0], pos[1], \
                .01 * numpy.cos(angle), .01 * numpy.sin(angle), \
                head_width=1, fc='none')

            self.background = None
            self.canvas.mpl_connect('draw_event', self.update_background)
            self.ax.set_xlim(-25, 25)
            self.ax.set_ylim(-25, 25)
            self.ax.set_xlabel('x')
            self.ax.set_ylabel('y')


            # open the figure window
            self.fig.canvas.draw()
            pl.show()


    def reset(self, pos=[0., 0., 0.], eul=[0., 0., 0.]):
        """Reset trajectory arrays."""
        print 'resetting'

        # store data
        self.posTraj = numpy.array(pos)
        self.eulTraj = numpy.array(eul)
        self.t = numpy.array(self.t_start)

        if self.realtimeplotting:

            if self.background is None: return False
            self.canvas.restore_region(self.background)

            self.line.set_data(self.posTraj[0], self.posTraj[1])
            self.ax.draw_artist(self.line)

            if self.headdirection:
                angle = numpy.pi / 2 - eul[0] * numpy.pi / 180
                pl.arrow(pos[0], pos[1], \
                    .01 * numpy.cos(angle), .01 * numpy.sin(angle), \
                    head_width=1, fc='none')

            self.canvas.blit(self.ax.bbox)

        return True

    def update_background(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)


    def update(self, pos=[0., 0., 0.], eul=[0., 0., 0.], t=-1):
        """Draw new data consisting of x and y."""
        # print 'updating'

        # store data
        self.posTraj = numpy.vstack((self.posTraj, pos))
        if eul:
            self.eulTraj = numpy.vstack((self.eulTraj, eul))

        if t >= 0.0:
            self.t = numpy.append(self.t, t)

        if self.realtimeplotting:

            if self.background is None: return False
            self.canvas.restore_region(self.background)

            self.line.set_data(self.posTraj[:, 0], self.posTraj[:, 1])
            self.ax.draw_artist(self.line)

            if self.headdirection:
                angle = numpy.pi / 2 - eul[0] * numpy.pi / 180
                pl.arrow(pos[0], pos[1], \
                    .01 * numpy.cos(angle), .01 * numpy.sin(angle), \
                    head_width=1, fc='none')

            self.canvas.blit(self.ax.bbox)

        return True

    def save(self, fileNamePart):
        if len(self.t.shape):
            data = numpy.hstack((numpy.vstack(self.t), self.posTraj))
        else:
            data = self.posTraj
        if len(data.shape) == 1:  # delete everything in case just the initialization value is in the array
            data = numpy.array([])
        else:
            data = data[1:]  # skip first entry, it was just for initialization
        save_trajectory(fileNamePart + '_position.traj', data, self.meta)

        if len(self.eulTraj):
            if len(self.t.shape):
                data = numpy.hstack((numpy.vstack(self.t), self.eulTraj))
            else:
                data = self.eulTraj
            if len(data.shape) == 1:  # delete everything in case just the initialization value is in the array
                data = numpy.array([])
            else:
                data = data[1:]  # skip first entry, it was just for initialization
            save_trajectory(fileNamePart + '_euler.traj', data, self.meta)


class taskTraj_recorder(traj_recorder):
    """Extended class to record task trajectories."""
    
    def __init__(self, pos=[0, 0, 0, 0, 0], eul=[0., 0., 0.], task='', meta={}, realtimeplotting=0):
        traj_recorder.__init__(self, pos, eul, meta, realtimeplotting)
        self.taskTraj = [task]

    def reset(self, pos=[0., 0., 0.], eul=[0., 0., 0.], task=''):
        traj_recorder.reset(self, pos, eul)
        self.taskTraj = [task]
    
    def update(self, pos=[0., 0., 0.], eul=[0., 0., 0.], task='', t=-1):
        traj_recorder.update(self, pos, eul, t)
        self.taskTraj.append(task)

    def save(self, fileNamePart):
        if len(self.t.shape):
            data = numpy.hstack((numpy.vstack(self.t), self.posTraj))
        else:
            data = self.posTraj
        data = numpy.hstack((data, numpy.vstack(self.taskTraj)))
        if len(data.shape) == 1:  # delete everything in case just the initialization value is in the array
            data = numpy.array([])
        else:
            data = data[1:]  # skip first entry, it was just for initialization
        save_trajectory(fileNamePart + '_position.traj', data, self.meta)

        if len(self.eulTraj):
            if len(self.t.shape):
                data = numpy.hstack((numpy.vstack(self.t), self.eulTraj))
            else:
                data = self.eulTraj
            data = numpy.hstack((data, numpy.vstack(self.taskTraj)))
            if len(data.shape) == 1:  # delete everything in case just the initialization value is in the array
                data = numpy.array([])
            else:
                data = data[1:]  # skip first entry, it was just for initialization
            save_trajectory(fileNamePart + '_euler.traj', data, self.meta)

    
class decisionTraj_recorder(traj_recorder):
    """Extended class to record trajectories related
        to binary decision tasks."""

    def __init__(self, pos=[0, 0, 0, 0, 0], eul=[0., 0., 0.], meta={}, realtimeplotting=0, realtimeSeqAna=0):
        traj_recorder.__init__(self, pos, eul, meta, realtimeplotting)

        if meta.has_key('task'):
            self.task = meta['task']
        if meta.has_key('taskStimulus'):
            self.taskStimulus = meta['taskStimulus']
        if meta.has_key('stimuli'):
            exec 'self.stimuli=' + meta['stimuli']

        self.realtimeSeqAna = realtimeSeqAna
        if self.realtimeSeqAna:
            self.__init_realtimeSeqAna()

    def __init_realtimeSeqAna(self):
        pl.ion()
        self.fig = pl.figure(facecolor='w')
        self.ax = self.fig.add_axes([.15, .12, .8, .8])
        # axSmall = self.fig.add_axes([.15, .85, .8, .07])
        self.canvas = self.ax.figure.canvas
        self.ax.grid()

        # draw test lines
        n = numpy.arange(50)
        alphas = [.05, .01, .001]
        h1_line = []
        h0_line = []
        for i, alpha in enumerate(alphas):
            h1_line.append(decisionLine(n, h=1, alpha=alpha))
            h0_line.append(decisionLine(n, h=0, alpha=alpha))

        for i, alpha in enumerate(alphas):
            self.ax.plot(n, h1_line[i], color=custom_plot.colors[i], linewidth=3, alpha=.75, label=str(alpha))
            self.ax.plot(n, h0_line[i], color=custom_plot.colors[i], linewidth=3, alpha=.75)
        # custom_plot.huebschMachen(self.ax)
        self.ax.set_xlabel('Decision #')
        self.ax.set_ylabel('Right decisions')
        self.ax.set_xlim(0, n[-1])
        self.ax.set_ylim(0, n[-1])
        self.ax.legend(loc='upper left')  # , frameon=False)

        # get line of data to be able to extend it later
        self.line = Line2D([0], [0], linestyle='-', color='b', marker='o')
        self.ax.add_line(self.line)

        self.background = None
        self.canvas.mpl_connect('draw_event', self.update_background)

        # small axis
        # if self.task=='visual':
        #    axSmall.axvspan(0.25, .75, color=self.stimuli[abs(self.taskStimulus-1)])
        #    axSmall.axvspan(1.25, 1.75, color=self.stimuli[self.taskStimulus])
        #    axSmall.text(0.5, 1.1, str(self.stimuli[abs(self.taskStimulus-1)]), horizontalalignment='center')
        #    axSmall.text(1.5, 1.1, str(self.stimuli[self.taskStimulus]), horizontalalignment='center')
        # axSmall.text(0.5, -.6, 'wrong', horizontalalignment='center')
        # axSmall.text(1.5, -.6, 'right', horizontalalignment='center')
        # custom_plot.allOff(axSmall)
        # axSmall.set_xlim(0, 2)

        # open the figure window
        self.fig.canvas.draw()
        pl.show()


    def getRightDecision(self):

        if self.posTraj[-1, 0]:
            chosenColor = self.posTraj[-1, 1]
        else:
            chosenColor = abs(self.posTraj[-1, 1] - 1)

        return int(self.taskStimulus == chosenColor)


    def realtimeSequentialAnalysis(self):
        pass


    def update_background(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)

    def update(self, pos=[0., 0., 0., 0., 0.], eul=[0., 0., 0.], t=-1):
        traj_recorder.update(self, pos, eul, t)

        if self.realtimeSeqAna:
            if self.background is None: return False
            self.canvas.restore_region(self.background)

            data = self.line.get_data()
            data0 = numpy.append(data[0], data[0][-1] + 1)
            data1 = numpy.append(data[1], data[1][-1] + self.getRightDecision())
            self.line.set_data(data0, data1)
            self.ax.draw_artist(self.line)
            self.canvas.blit(self.ax.bbox)

    def save(self, fileNamePart):
        traj_recorder.save(self, fileNamePart)
        self.fig.savefig(fileNamePart + 'SeqAna.pdf', format='pdf')
