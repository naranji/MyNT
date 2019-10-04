"""
A module for spatial paths, i.e., trajectories.
"""

__author__ = ("KT", "Moritz Dittmeyer", "Benedikt Ludwig", "Sven Schoernich")
__version__ = "4.61, July 2012"

# python modules
import sys, struct, os

# other modules
import numpy
import scipy.stats

import matplotlib.pyplot as pl
from matplotlib.collections import LineCollection
import matplotlib as mpl

# custom made modules
import custom_plot, signale


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

    # in case there are no values in the file
    if traj.size < 5 and traj.mean() == 0:
        traj = numpy.array([[0., 0., 0., 0.]])

    if meta['mazetype'] == 'linearMaze':
        return VRlinearMazeTrajectory(traj, meta)
    elif meta['mazetype'].find('yDecision') + 1:
        return decisionMazeTrajectory(traj, meta)
    else:
        return VRtrajectory(traj, meta)


def load_collisionTrajectory(fileName):

    traj = numpy.loadtxt(fileName)
    if traj.shape.__len__() == 1:
        traj = traj.reshape(1, traj.shape[-1])

    # in case there are no values in the file
    if traj.size < 5 and traj.mean() == 0:
        traj = numpy.array([[0., 0., 0., 0.]])

    meta = _read_metadata(fileName)
    if meta['mazetype'] == 'linearMaze':
        return linearMazeCollisionTrajectory(traj, meta)
    else:
        return collisionTrajectory(traj, meta)



def load_nvtFile(fileName, mazeType='', showHeader=False):
    """
    For loading Neuralynx video tracker data
    with python/numpy. Returns a trajectory object.
    """
    timeStamps = []
    posSamples = []
    eulSamples = []

    with open(fileName, 'rb') as f:

        header = f.read(16 * 2 ** 10)  # the first 16 kB are an ASCII text header

        if showHeader:
            print header

        count = 0
        while True:
            swstx = f.read(2)  # UInt16, Value indicating the beginning of a record.
                                        # Always 0x800 (2048).
            swid = f.read(2)  # UInt16, ID for the originating system of this record.
            swdata_size = f.read(2)  # UInt16, Size of a VideoRec in bytes.
            qwTimeStamp = f.read(8)  # UInt64, Cheetah timestamp for this record.
                                        # This value is in microseconds.

            dwPoints = f.read(400 * 4)  # UInt32[], Points with the color bitfield values
                                        # for this record. This is a 400 element array.
                                        # See Video Tracker Bitfield Information.

            sncrc = f.read(2)  # Int16, Unused*
            dnextracted_x = f.read(4)  # Int32, Extracted X location of the target being tracked.
            dnextracted_y = f.read(4)  # Int32, Extracted Y location of the target being tracked.
            dnextracted_angle = f.read(4)  # Int32, The calculated head angle in degrees clockwise
                                            # from the positive Y axis.
                                            # Zero will be assigned if angle tracking is
                                            # disabled.
            dntargets = f.read(50 * 4)  # Int32[], Colored targets using the same bitfield format
                                        # used by the dwPoints array.
                                        # Instead of transitions, the bitfield indicates
                                        # the colors that make up each particular target and
                                        # the center point of that target.
                                        # This is a 50 element array sorted by
                                        # size from largest (index 0) to smallest (index 49).
                                        # A target value of 0 means that no target is present in
                                        # that index location.
                                        # See Video Tracker Bitfield Information below.

            if qwTimeStamp:
                qwTimeStamp = struct.unpack('L', qwTimeStamp)[0]
                dnextracted_x = struct.unpack('i', dnextracted_x)[0]
                dnextracted_y = struct.unpack('i', dnextracted_y)[0]
                dnextracted_angle = struct.unpack('i', dnextracted_angle)[0]

                dwPoints = [dwPoints[i:i + 4] for i in range(0, dwPoints.__len__(), 4)]
                points = []
                for point in dwPoints:
                    point = struct.unpack('I', point)[0]
                    points.append(point)
                dwPoints = numpy.array(points)

                dntargets = [dntargets[i:i + 4] for i in range(0, dntargets.__len__(), 4)]
                targets = []
                for target in dntargets:
                    target = struct.unpack('i', target)[0]
                    targets.append(target)
                dntargets = numpy.array(targets)

                timeStamps.append(qwTimeStamp)
                posSamples.append([dnextracted_x, dnextracted_y, 0])
                eulSamples.append([dnextracted_angle, 0, 0])


                # print some data?
                if count < 0:
                    count += 1
                    print '----> count', count
                    print qwTimeStamp
                    print dnextracted_x
                    print dnextracted_y
                    print dnextracted_angle
                    # print dwPoints
                    print dwPoints.shape
                    # print dntargets
                    print dntargets.shape
                    print ''
            else:
                break

    timeStamps = numpy.array(timeStamps) / 1000. / 1000.  # change to s,
                                                        # Neuralynx time stamps are in us
    posSamples = numpy.array(numpy.hstack((numpy.vstack(timeStamps), posSamples)))
    eulSamples = numpy.array(numpy.hstack((numpy.vstack(timeStamps), eulSamples)))



    filePath = os.path.dirname(os.path.abspath(fileName))
    meta = {'file' : fileName, 'path' : filePath}
    if 'VT_parameters.dat' in os.listdir(filePath):  # load parameters.dat
        meta.update(signale._read_metadata(filePath + '/VT_parameters.dat'))

    if mazeType == 'linearMaze':
        meta.update({'mazetype' : mazeType})
        posTraj = linearMazeTrajectory(posSamples, meta)
        eulTraj = []
        # eulTraj = linearMazeTrajectory(eulSamples, meta)
    else:
        posTraj = trajectory(posSamples, meta)
        eulTraj = trajectory(eulSamples, meta)

    return posTraj, eulTraj




def load_rewardTrajectory(fileName):

    traj = numpy.loadtxt(fileName)
    if traj.shape.__len__() == 1:
        traj = traj.reshape(1, traj.shape[-1])
    meta = _read_metadata(fileName)

    # in case there are no values in the file
    if traj.size < 6 and traj.mean() == 0:
        traj = numpy.array([[0., 0., 0., 0., 0.]])

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


def load_stimuli(stimuliDates, folderNames, rewardTraj, folderName):
    """For loading the cvs stimuli files for the yDecisionwCues experiments.

        Rather specifically tuned to the yDecisionwCues experiments.
        I don't really like it to be in here,butalso don't have a better place at the time and
        I need to use it several times.
    """
    stimuli = []
    stimuliInFile = []
    folderNames_dummy = numpy.round(folderNames / 100) - 20000000  # change folderNames array to cope with stimuli date style
    for i, d in enumerate(stimuliDates):
        stimuliFile = [f for f in  os.listdir(folderName) if f.find(d) + 1][0]  # get stimulus file

        # load stimuli file
        if stimuliFile.endswith('.csv'):
            stimuliInFile.append(numpy.loadtxt(folderName + stimuliFile, delimiter=',', dtype='str'))  # csv fromat => comma delimiters
        elif stimuliFile.endswith('.tsv'):
            stimuliInFile.append(numpy.loadtxt(folderName + stimuliFile, delimiter='\t', dtype='str'))  # csv fromat => comma delimiters

        # get number of decisions made with corresponding stimuliFile
        if i + 1 < len(stimuliDates):
            indices1, = numpy.where(folderNames_dummy < int(stimuliDates[i + 1]))
            indices2, = numpy.where(folderNames_dummy >= int(stimuliDates[i]))
            indices = numpy.intersect1d(indices1, indices2)
        else:
            indices, = numpy.where(folderNames_dummy >= int(d))
        numDecisions = 0
        for index in indices:
            numDecisions += rewardTraj[index].numPlaces  # numPlaces in a decisions trajectory is the number of decisions
        print 'more decisions than stimuli?', numDecisions > stimuliInFile[-1].shape[0], numDecisions, stimuliInFile[-1].shape[0]
        stimuliInFile[-1] = stimuliInFile[-1][:numDecisions]  # only include stimuli with appropriate number
        stimuli.extend(stimuliInFile[-1].tolist())
    stimuli = numpy.array(stimuli)

    # get names of the stimulus sets
    stimulusSetsNames = stimuli[:, 0]
    stimulusSetsNames = dict.fromkeys(stimulusSetsNames).keys()
    stimulusSetsNames.sort()
    # would be finished here, if there are no negative 'names'
    # but for negative names -1 is lower than -2
    if stimuliFile.endswith('.csv'):
        stimulusSetsNames = sorted([int(s[:-1]) for s in stimulusSetsNames[::2]])
        dummy = []
        for s in stimulusSetsNames:
            dummy.append(str(s) + 'a')
            dummy.append(str(s) + 'b')
        stimulusSetsNames = dummy

    # get the different stimulus sets
    stimulusSets = []
    for s in stimulusSetsNames:
        stimulusSets.append(stimuli[numpy.where(stimuli[:, 0] == s)[0][0]])
    stimulusSets = numpy.array(stimulusSets)



    # analyze stimulus files
    for i, f in enumerate(stimuliInFile):
        decisions_sum_perSubSet = numpy.array([])
        for set in stimulusSetsNames:
            indices = numpy.where(f[:, 0] == set)[0]
            decisions_sum_perSubSet = numpy.append(decisions_sum_perSubSet, indices.shape[0])
        print 'file:', stimuliDates[i]
        # print 'presentations of sets 0-3 (%):', numpy.round(numpy.sum(decisions_sum_perSubSet[:-2])/f.shape[0], 2)
        # print 'presentations of set 4 (%):', numpy.round(numpy.sum(decisions_sum_perSubSet[8:])/f.shape[0], 2)
        print 'number of decisions in file:', f.shape[0]


    return stimuli, stimulusSetsNames, stimulusSets


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



###################################################### CLASSES


class trajectory(object):

    def __init__(self, traj, meta):

        # initialize meta data
        self.meta = meta
        self.dt = 0.1  # claim a dt [s]
        self.t_start = 0.0
        self.mazeType = ''
        self.timeUnit = 's'
        self.spaceUnit = 'm'
        self.rewardsArea = []
        self.view_eulerOffsets = [0, 0, 0]
        self.euler = -numpy.array(self.view_eulerOffsets)  # for storing the current euler angles
        self.euler %= 360


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
        print "dt is", self.dt


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
            laps, lapTimes = self.getLaps()
            for item in lapTimes:
                lapStart = self.getIndexFromTime(item[0])
                lapEnd = self.getIndexFromTime(item[1])
                diffsTmp.append(numpy.diff(self.places[lapStart:lapEnd], axis=0)[:, 0:2])

            diffs = []
            for item in diffsTmp:
                for subitem in item:
                    diffs.append(subitem)
            diffs = numpy.array(diffs)

        else:
            diffs = numpy.diff(self.places, axis=0)[:, 0:2]


        if vec:
            speed = diffs / self.dt
        else:
            speed = numpy.sqrt(diffs[:, 1] ** 2 + diffs[:, 0] ** 2) / self.dt  # [space units/s]

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
        acceleration = diffs / self.dt  # [space units/s^2]

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
        jerk = diffs / self.dt  # [space units/s^3]

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
            index = self.getIndexFromTime(time)  # get index of the time point that is just smaller than time
            place.extend(self.places[[index]])
            place = numpy.array(place)

        return place

    def getRunningTraj(self, threshspeed=0.01, window_len=51):
        """ Get times & places, where the animal was running.

        threshspeed ... minimum running speed
        """

        # get running speed and smooth it a bit
        speed_dummy = signale.smooth(self.getSpeed()[1], window_len=window_len, window='hanning')

        indices = numpy.where(speed_dummy < threshspeed)[0]
        self.run_times = numpy.delete(self.times, indices, 0)
        self.run_places = numpy.delete(self.places, indices, 0)

        return self.run_times, self.run_places


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

        if not fig:
            fig = pl.figure(figsize=(8, 6))
        if not ax:
            ax = fig.add_subplot(111)

        line_segments = LineCollection([[x, self.places[i + 1 + offset, 0:2]] \
                            for i, x in enumerate(self.places[offset:-(1 + offset), 0:2])], \
                            linestyles='solid', linewidths=mpl.rcParams['lines.linewidth'] / 2.)
        line_segments.set_array(self.times)
        ax.add_collection(line_segments)

        axcb = fig.colorbar(line_segments)
        ax.plot(self.places[0 + offset, 0], self.places[0 + offset, 1], 'o')  # start point
        ax.plot(self.places[-2, 0], self.places[-2, 1], 'd')  # end point


        # huebsch machen
        custom_plot.huebschMachen(ax)
        if chic:
            pass
        else:
            ax.set_title(self.mazeType)
        if language == 'd':  # in german
            axcb.set_label('Zeit (' + self.timeUnit + ')')
            ax.set_xlabel('x-Position (' + self.spaceUnit + ')')
            ax.set_ylabel('y-Position (' + self.spaceUnit + ')')
        else:  # in english
            axcb.set_label('Time (' + self.timeUnit + ')')
            ax.set_xlabel('x position (' + self.spaceUnit + ')')
            ax.set_ylabel('y position (' + self.spaceUnit + ')')
        dx = numpy.round((int(self.t_stop) / 60 * 60 - int(self.t_start) / 60 * 60) / 4.)
        if not dx:
            dx = 60
        xticks = numpy.arange(round(self.t_start), round(self.t_stop) + 1, dx)
        axcb.set_ticks(xticks)

        xoffset = self.xWidth * .15
        yoffset = self.yWidth * .15

# #        xmin = -xoffset
# #        xmax = self.xWidth + xoffset
# #        ymin = -yoffset
# #        ymax = self.yWidth + yoffset

        xmin = self.xlim[0] - xoffset
        xmax = self.xlim[1] + xoffset
        ymin = self.ylim[0] - yoffset
        ymax = self.ylim[1] + yoffset

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


        if chic:
            custom_plot.allOff(ax)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')
# #            pos1[3] = .45
# #            pos2[3] = .45
# #            ax1.set_position(pos1)
# #            ax2.set_position(pos2)

            minX = self.xlim[0]
            maxX = self.xlim[1]
            dX = numpy.round((self.xWidth) / 10, 1)

            minY = self.ylim[0]
            maxY = self.ylim[1]
            dY = numpy.round((self.yWidth) / 10, 1)

            ax.hlines(y=minY - dY / 2 , xmin=maxX - dX , xmax=maxX)
            ax.vlines(x=minX - dX / 4 , ymin=minY , ymax=minY + dY)

            print "horizontal scale bar: ", dX
            print "vertical scale bar: ", dY



        if not fig:
            pl.show()

        return fig, ax, axcb


    def plotCumPlaces(self, fig=None, ax=None):

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        cumPlace = self.cumPlaces()
        ax.plot(self.times, cumPlace)

        ax.set_xlabel('Time (' + self.timeUnit + ')')
        ax.set_ylabel('Cumulative path length (' + self.spaceUnit + ')')
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
                index_min = find_nearest_vec(self.places[index2:index1], places[i])[0] + index2
                time_diffs = numpy.append(time_diffs, t - self.times[index_min])
            self.time_offset(time_diffs.mean())
            # print "time differences:", time_diffs
            print "=> shifted time axis by" , time_diffs.mean(), ", std is" , time_diffs.std()
            self.__recalc_startstop()


    def time_offset(self, offset=0.):
        self.times += offset
        self.__recalc_startstop()


class VRtrajectory(trajectory):

    def __init__(self, traj, meta):
        trajectory.__init__(self, traj, meta)

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

        if meta.has_key('yaw'):
            self.yaw = meta['yaw'] % 360
            self.euler[0] += self.yaw
            self.euler %= 360

        self.oriented = False
        self.orient()

    def centerTraj(self):
        """
        Shifts the trajectory center to position (0,0).
        """
        trajectory.centerTraj(self)

        for i, p in enumerate(self.rewardsPos):
            self.rewardsPos[i, 0] += self.xWidth / 2
            self.rewardsPos[i, 1] += self.yWidth / 2


class linearMazeTrajectory(trajectory):

    def __init__(self, places, meta, initialized=False):
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
            places = self.run_places
            times = self.run_times
        else:
            places = self.places
            times = self.times

        xdir = numpy.diff(signale.smooth(places[:, 0], window_len=11, window='hanning'), axis=0)  # x-richtung bestimmen!

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


# #    def orient(self):
# #        """
# #        Orients the trajectory to yaw = 0 deg by projecting it
# #        via a dot product with the slope vector.
# #        """
# #        if not self.oriented:
# #            self.turn(-self.euler[0])
# #
# #            # rotate to positive axes, if necessary
# #            if self.yaw > 90 and self.yaw <= 270:
# #                self.places[:, 0:2] *= -1
# #
# #            self.getTrajDimensions()
# #            self.trackWidth = self.yWidth
# #            self.oriented = True
# #        else:
# #            print "NOTE: Track already oriented."


    def plot(self, fig=None, ax=None, offset=2, language='e', chic=False):

        fig, ax, axcb = trajectory.plot(self, fig, ax, offset, language, False)

        # adjust axes
        pos = [.15, .63, .65, .3]
        ax.set_position(pos)
        poscb = list(axcb.ax.get_position().bounds)
        poscb[0] = .825
        poscb[1] = pos[1]
        poscb[3] = pos[3]
        axcb.ax.set_position(poscb)


        # huebsch machen
        xoffset = self.xWidth * .1
        yoffset = self.yWidth

# #        xmin = -xoffset
# #        xmax = self.xWidth + xoffset
# #        ymin = -yoffset
# #        ymax = self.yWidth + yoffset

        xmin = self.xlim[0] - xoffset
        xmax = self.xlim[1] + xoffset
        ymin = self.ylim[0] - yoffset
        ymax = self.ylim[1] + yoffset

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)


        # put additional axes
        pos1 = list(pos)
        pos1[1] = .15
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

            minX = self.xlim[0]
            maxX = self.xlim[1]
            dX = numpy.round((self.xWidth) / 10, 1)

            minY = self.ylim[0]
            maxY = self.ylim[1]
            dY = numpy.round((self.yWidth) / 2, 1)

            ax.hlines(y=minY - dY / 2 , xmin=maxX - dX , xmax=maxX)
            ax.vlines(x=minX - dX / 4 , ymin=minY , ymax=minY + dY)

            print "horizontal scale bar: ", dX
            print "vertical scale bar: ", dY


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

    def plotSpeedvsComponents(self, thresh=None, fig=None, ax=None, color='k', lineStyle='.'):

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        time, speed = self.getSpeed()  # to do: including thresh
        x = self.getXComponents()[1:]
        y = self.getYComponents()[1:]


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

        return fig, ax


class VRlinearMazeTrajectory(VRtrajectory, linearMazeTrajectory):

    def __init__(self, places, meta, initialized=False):
        if not initialized:
            VRtrajectory.__init__(self, places, meta)

        if not meta.has_key('expType'):
            meta['expType'] = 'VR'
        linearMazeTrajectory.__init__(self, places, meta, initialized=True)

        # get inter reward area distance
        # the distance between the reward positions.
        # The reward positions are expected to be at the track ends.
        self.IRAreaDistance = numpy.diff(self.rewardsPos, axis=0)
        self.IRAreaDistance = numpy.sqrt(self.IRAreaDistance[0, 0] ** 2 + self.IRAreaDistance[0, 1] ** 2)
        self.IRAreaDistance -= 2 * self.rewardsArea[0]
        # NOTE: previous command assumes circular reward zone with radius rewardsArea[0]


    def orient(self):
        """
        Orients the trajectory to yaw = 0 deg by projecting it
        via a dot product with the slope vector.
        """
        linearMazeTrajectory.orient(self)

        # rotate to rewardsPos positive axes, if necessary
        if self.yaw > 90 and self.yaw <= 270:
            self.rewardsPos[:, 0:2] *= -1

        # get inter reward area distance
        # the distance between the reward positions.
        # The reward positions are expected to be at the track ends.
        self.IRAreaDistance = numpy.diff(self.rewardsPos, axis=0)
        self.IRAreaDistance = numpy.sqrt(self.IRAreaDistance[0, 0] ** 2 + self.IRAreaDistance[0, 1] ** 2)
        self.IRAreaDistance -= 2 * self.rewardsArea[0]
        # NOTE: previous command assumes circular reward zone with radius rewardsArea[0]


    def plot(self, fig=None, ax=None, offset=2, language='e', chic=False):

        fig, ax, axcb = linearMazeTrajectory.plot(self, fig, ax, offset, language, chic)

        # draw rewardsPos
        for rPos in self.rewardsPos:
            if self.rewardsAreaType == 'circular':
                circ = pl.Circle(rPos, radius=self.rewardsArea[0], color=custom_plot.grau2)
                ax.add_artist(circ)
                custom_plot.drop_shadow_patches(ax, circ, sigma=5, offsets=(2, 2))
            elif self.rewardsAreaType == 'rectangular':
                pass

        # huebsch machen
        xoffset = self.xWidth * .1
        yoffset = self.yWidth

# #        xmin = -xoffset
# #        xmax = self.xWidth + xoffset
# #        ymin = -yoffset
# #        ymax = self.yWidth + yoffset

        xmin = self.xlim[0] - xoffset
        xmax = self.xlim[1] + xoffset
        ymin = self.ylim[0] - yoffset
        ymax = self.ylim[1] + yoffset

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        return fig, ax, axcb



class decisionMazeTrajectory(VRtrajectory):

    def __init__(self, places, meta):
        VRtrajectory.__init__(self, places, meta)

        self.wait_times = []
        if meta.has_key('wait_time'):
            self.wait_times.append(meta['wait_time'])
        if meta.has_key('wait_time_wrong'):
            self.wait_times.append(meta['wait_time_wrong'])
        self.wait_times = numpy.array(self.wait_times)


    def plot(self, times=[], fig=None, ax=None, offset=2, language='e', chic=False):

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


class collisionTrajectory(VRtrajectory):

    def __init__(self, traj, meta):
        VRtrajectory.__init__(self, traj, meta)
        # self.times=traj[:,0]

    def getICI(self, ICImin=-1, t_start=0., t_stop=None):  # get Intercollision Interval
        """ Inter collision interval.
        """
        if ICImin == -1:
            ICImin = self.dt

        times = numpy.append(t_start, self.times)
        if t_stop:
            times = numpy.append(times, t_stop)

        ICI = numpy.diff(times)
        self.reducedICI = numpy.array([])
        # self.reducedPlaces = numpy.array([])
        # self.reducedTimes = numpy.array([])
        for i, ic in enumerate(ICI):
            if ic >= ICImin:
                self.reducedICI = numpy.append(self.reducedICI, ic)
                # self.reducedPlaces = numpy.append(self.reducedPlaces, self.places[i+1])
                # self.reducedTimes = numpy.append(self.reducedTimes, times[i+1])
        return self.reducedICI

    def getICD(self, traj, ICDmin=0, t_start=0., t_stop=None, withTime=False):  # get Intercollision distance
        """Inter collision distance.
        """
        indices = numpy.array([])
        distances = numpy.array([])
        distances_time = numpy.array([])
        cumDistances = numpy.array([])

        times = numpy.append(t_start, self.times)
        if t_stop:
            times = numpy.append(times, t_stop)

        for time in times:
            index = traj.getIndexFromTime(time)
            indices = numpy.append(indices, index)
        for i, index in enumerate(indices[0:-1]):
            if index < indices[i + 1]:  # just a sanity check
                d = traj.cumPlaces(index, indices[i + 1])
                if d[-1] >= ICDmin:
                    cumDistances = numpy.append(cumDistances, d)
                    distances = numpy.append(distances, d[-1])
                    if withTime:
                        distances_time = numpy.append(distances_time, times[i + 1])
        if withTime:
            return distances_time, distances, cumDistances
        else:
            return distances, cumDistances

    def plotICI(self, ICImin=-1, fig=None, ax=None):

        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        ICI = self.getICI(ICImin)
        ax.plot(ICI)

        ax.set_xlabel('Collision #')
        ax.set_ylabel('Intercollision interval [s]')
        pl.show()

        return fig, ax

    def plotICD(self, traj, ICDmin=0, fig=None, ax=None):
        [ICD_times, ICD, cumICD] = self.getICD(traj, ICDmin, withTime=True)

        if not fig:
            fig = pl.figure()
        # if not ax:
        #    ax = fig.add_subplot(111)
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

        return fig, ax

    def removeCollisions(self, center=[0, 0], radius=0):
        """Remove collisions inside the cicular area of radius around center."""

        diffs = self.places[:, 0:2] - center
        distance = numpy.sqrt(diffs[:, 1] ** 2 + diffs[:, 0] ** 2)
        indices = numpy.where(distance < radius)
        if indices[0].size:
            self.times = numpy.delete(self.times, indices)
            self.places = numpy.delete(self.places, indices, 0)
        self.numPlaces = self.places.shape[0]

    def removeCollisionsAtRewardPositions(self):
        """Remove collisions at reward positions."""
        for r in self.rewardsPos:
            self.removeCollisions(r, self.rewardsArea[0])

    def purgeCollisions(self, dt=None):
        """Delete double entries (i.e., too close entries from collision trajectory)."""
        if not dt:
            dt = self.dt

        indices = numpy.where(numpy.diff(self.times) < dt)
        if indices[0].size:
            self.times = numpy.delete(self.times, indices)
            self.places = numpy.delete(self.places, indices, 0)
        self.numPlaces = self.places.shape[0]


class linearMazeCollisionTrajectory(collisionTrajectory, VRlinearMazeTrajectory):

    def __init__(self, traj, meta):
        collisionTrajectory.__init__(self, traj, meta)
        VRlinearMazeTrajectory.__init__(self, self.places, meta, initialized=True)
        # self.times=traj[:,0]

    def removeCollisionsAtRewardPositions(self):
        """Remove collisions at reward positions."""

        before = self.times.size
        collisionTrajectory.removeCollisionsAtRewardPositions(self)
        print self.times.size, "collisions along track ("\
            + str(100 * self.times.size / before), "%) from", \
            before, "including reward sites"


    def collisionSorting(self, traj, alphamax=15.):
        """
        alphamax ... deg, maximum angle relative to wall not to be counted as collision
        """

        # definitions and initializations
        alphamaze = traj.euler[0]  # orientation of the maze in Vizard's left-handed coordinate system

        collisionTimes = []
        collisionPlaces = []

        offset = 3  # index offset of position before collision
        nanOffset = 5  # index offset is alpha is declared NaN (i.e. confirm if really NaN)


# #        # transform to right-handed coordinate system:
# #        if 0. <= alphamaze <= 180.:
# #            alphamaze = abs(180. - alphamaze)
# #        elif 180. < alphamaze < 360.:
# #            alphamaze = 360. - abs(180. - alphamaze)
# #        else:
# #            print "Maze orientation yaw larger than 360 deg! \nCalculation still to be implemented!"


        # keep orientation between 0 and 180 degrees:
        if alphamaze > 180.:
            alphamaze = alphamaze - 180.
        # print "yaw of the maze", alphamaze

        # calculate collision angles and check if collision fulfills criteria
        for i, item in enumerate(self.times):

            k = traj.getIndexFromTime(item)

            # calculates (probably) correct index of trajectory at collision
            getCorrectIndex = True
            if getCorrectIndex:
                # find correct collision index
                # Idea: if translation is very low between x0 and x1, there probably
                # was a collision at x0

                indices = []
                lower = k - 1
                upper = k + 1

                # exceptions if index out of range
                if k == 0:
                    lower = k

                kmax = traj.getIndexFromTime(self.times[-1])
                if k == kmax:
                    upper = k - 1
                elif k == kmax - 1:
                    upper = k

                for j in range(lower , upper + 1):
                    indices.append([j , numpy.linalg.norm(traj.places[j + 1] - traj.places[j])])
                index = indices[numpy.argmin(indices, 0)[1]][0]
                k = numpy.int(index)


            xpos = []
            ypos = []

            if k == traj.getIndexFromTime(self.times[-1]):
                k_uplim = k
            else:
                k_uplim = k + 1  # in general, take the position one index after calculated collision into account
                                # Is that a good choice after using KillCollisionBuffer?

            for j in range(k - offset, k_uplim + 1):
                xpos.append(traj.places[j, 0])
                ypos.append(traj.places[j, 1])
                j += 1

            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xpos, ypos)
            alpha = numpy.arctan(slope) * 180. / numpy.pi



            # if angle is NaN, confirm if true by using a larger offset
            if numpy.isnan(alpha):

                j = k - nanOffset
                ypos = []
                xpos = []
                for j in range(k - nanOffset , k_uplim + 1):
                    xpos.append(traj.places[j, 0])
                    ypos.append(traj.places[j, 1])
                    j += 1

                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xpos, ypos)
                alpha = numpy.arctan(slope) * 180. / numpy.pi


            # make sure the relative angle is processed correctly
            alpharel = abs(alphamaze - alpha)


            # use only sharp angles
            if alpharel > 270.:
                alpharel = abs(360. - alpharel)
            elif alpharel > 90.:
                alpharel = abs(180. - alpharel)


            # check if collision is true (i.e. angle larger than maximum non-collision angle)
            if alpharel > alphamax and alpharel != 0.:
                collisionTimes.append(self.times[i])
                collisionPlaces.append(self.places[i])
                # print "index", i, "alpha", alpha, "alpharel", alpharel


        print collisionTimes.__len__(), "true collisions (angle <=", \
                alphamax, "deg) remaining from", self.times.__len__()
        self.times = numpy.array(collisionTimes)
        self.places = numpy.array(collisionPlaces)



class rewardTrajectory(VRtrajectory):

    def __init__(self, traj, meta):
        if traj.shape[1] == 2:  # in case only times and IDs were stored fill up
            traj = numpy.hstack((traj, numpy.ones([traj.shape[0], 3])))
            print "NOTE: Only times and IDs in reward trajetory provided."
        dummy = numpy.hstack((traj[:, 0:1], traj[:, 2:]))
        VRtrajectory.__init__(self, dummy, meta)
        self.IDs = traj[:, 1]

    def getIRI(self, t_start=0., t_stop=None):  # get Interreward interval
        "Inter reward interval"
        times = numpy.append(t_start, self.times)
        if t_stop:
            times = numpy.append(times, t_stop)
        return numpy.diff(times)

    def getIRD(self, traj, t_start=0., t_stop=None):  # get Interreward distance
        """Get the inter reward distance

        Parameters:
        traj ... the corresponding normal trajectory

        Returns:
        distances, cumDistances

        """
        indices = numpy.array([])
        distances = numpy.array([])
        cumDistances = numpy.array([])

        times = numpy.append(t_start, self.times)
        if t_stop:
            times = numpy.append(times, t_stop)

        for time in times:
            index = traj.getIndexFromTime(time)
            indices = numpy.append(indices, index)
        for i, index in enumerate(indices[0:-1]):
            d = traj.cumPlaces(index, indices[i + 1])
            cumDistances = numpy.append(cumDistances, d)
            if d.size:
                distances = numpy.append(distances, d[-1])
        return distances, cumDistances

    def plotIRI(self, fig=None, ax=None):
        if not fig:
            fig = pl.figure()
        if not ax:
            ax = fig.add_subplot(111)

        IRI = self.getIRI()
        ax.plot(IRI)

        ax.set_xlabel('Reward #')
        ax.set_ylabel('Inter reward interval [s]')
        pl.show()

        return fig, ax


    def plotIRD(self, traj, fig=None, ax=None):
        """Plot the inter reward distance

        Parameters:
        traj ... the corresponding normal trajectory

        """
        [IRD, cumIRD] = self.getIRD(traj)

        if not fig:
            fig = pl.figure()
        # if not ax:
        #    ax = fig.add_subplot(111)
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        ax.plot(IRD)
        ax2.plot(cumIRD)

        # draw track length
        if self.mazeType == 'linearMaze':
            ax.plot([0, IRD.shape[0]], [self.IRAreaDistance, self.IRAreaDistance])
            ax2.plot([0, cumIRD.shape[0]], [self.IRAreaDistance, self.IRAreaDistance])

        ax2.set_xlabel('Reward #')
        ax.set_ylabel('Inter reward distance (' + self.spaceUnit + ')')
        ax2.set_ylabel('Cummulated distance btw. rewards (' + self.spaceUnit + ')')
        pl.show()

        return fig, ax


    def purgeRewards(self, dt=None):
        """Delete double entries (i.e., too close entries from reward trajectory)."""
        if not dt:
            dt = self.dt

        indices = numpy.where(numpy.diff(self.times) < dt)
        if indices[0].size:
            self.times = numpy.delete(self.times, indices)
            self.places = numpy.delete(self.places, indices, 0)

        self.numPlaces = self.places.shape[0]


class linearMazeRewardTrajectory(rewardTrajectory, VRlinearMazeTrajectory):

    def __init__(self, traj, meta):
        rewardTrajectory.__init__(self, traj, meta)
        VRlinearMazeTrajectory.__init__(self, self.places, meta, initialized=True)
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
        # the third column is the state (binary) of the left side
        self.leftSide_StimulusIndex = traj[:, 2]
        self.rightSide_StimulusIndex = abs(traj[:, 2] - 1)

    def getRightDecisions(self, display=False):
        """ Returns an array with the outcome of each decision.
        1 for right decision. 0 for wrong decision."""

        rightDecisions = []
        for i, id in enumerate(self.IDs):
            if id:
                chosenColor = self.leftSide_StimulusIndex[i]
            else:
                chosenColor = self.rightSide_StimulusIndex[i]
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
            patches = ax.bar(bins, hist / hist.sum() * 100, width=bar_width, \
                            color=custom_plot.grau, edgecolor=custom_plot.grau)
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

