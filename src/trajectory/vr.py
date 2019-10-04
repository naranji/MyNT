"""
trajectory.real
===============

A module for trajectories from virtual reality.
"""

__author__ = ("KT", "Benedikt Ludwig", "Sven Schoernich")
__version__ = "6.2, July 2013"


# python modules
import sys

# other modules
import numpy
import scipy.stats
import matplotlib.pyplot as pl
from matplotlib.collections import LineCollection

# custom made modules
import custom_plot

# package modules
import trajectory
from trajectory import trajectory
from trajectory import linearMazeTrajectory


###################################################### FUNCTIONS

###################################################### CLASSES

class VRtrajectory(trajectory):


    def __init__(self, traj, meta={}):
        trajectory.__init__(self, traj, meta)


        self.rewardsPos = []
        self.rewardsArea = []

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

        if meta.has_key('view_eulerOffsets'):
            exec 'self.view_eulerOffsets=' + meta['view_eulerOffsets']
            self.euler = -numpy.array(self.view_eulerOffsets)  # for storing the current euler angles
            self.euler %= 360
            self.turn(self.view_eulerOffsets[0])  # remove the yaw component of view_eulerOffsets



    def centerTraj(self):
        """
        Shifts the trajectory center to position (0,0).
        """
        trajectory.centerTraj(self)

        for i, p in enumerate(self.rewardsPos):
            self.rewardsPos[i, 0] += self.xWidth / 2
            self.rewardsPos[i, 1] += self.yWidth / 2



class VRlinearMazeTrajectory(VRtrajectory, linearMazeTrajectory):

    def __init__(self, places, meta={}, initialized=False):
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


    def getLaps(self):
        """
        Returns the number of laps in the trajectory.

        A lap is defined as the path/time between leaving a reward
        area and entering the next. NOTE, that the next reward area
        might be the same or a different reward area.
        """
        x = self.getXComponents()
        i = 0
        laps = 0
        inLap = False
        lapTimes = []
        while i < x.shape[0]:
            if x[i] > self.rewardsPos[0, 0] + self.rewardsArea[0] and x[i] < self.rewardsPos[1, 0] - self.rewardsArea[0] and not inLap:
                laps += 1
                lapStart = self.times[i]
                inLap = True
            elif inLap and (x[i] < self.rewardsPos[0, 0] + self.rewardsArea[0] or x[i] > self.rewardsPos[1, 0] - self.rewardsArea[0]):
                lapTimes.append([lapStart, self.times[i - 1]])
                inLap = False
            i += 1
        return laps, lapTimes


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

        xmin = -xoffset
        xmax = self.xWidth + xoffset
        ymin = -yoffset
        ymax = self.yWidth + yoffset

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        return fig, ax, axcb


class decisionMazeTrajectory(VRtrajectory):

    def __init__(self, places, meta={}):
        VRtrajectory.__init__(self, places, meta)
        if meta.has_key('yaw'):
            self.yaw = meta['yaw'] % 360
            self.euler[0] += self.yaw
            self.euler %= 360

        self.oriented = False

        self.wait_times = []
        if meta.has_key('wait_time'):
            self.wait_times.append(meta['wait_time'])
        if meta.has_key('wait_time_wrong'):
            self.wait_times.append(meta['wait_time_wrong'])
        self.wait_times = numpy.array(self.wait_times)


    def getLaps(self, times=[], offset=2):
        """
        Returns the number of laps in the trajectory.

        A lap is defined as the path/time between leaving a reward
        area and entering the next. NOTE, that the next reward area
        might be the same or a different reward area.
        """

        if self.wait_times.size:
            wait_time = self.wait_times[0]
        else:
            wait_time = 2.

        indices1 = numpy.array([offset])
        indices2 = numpy.array([])
        indices1 = numpy.hstack((indices1, self.getIndexFromTime(times + wait_time) + 1))
        indices2 = numpy.hstack((indices2, self.getIndexFromTime(times), self.times.size - 1))
        for i, idx1 in enumerate(indices1[:-1]):
            if indices1[i + 1] <= idx1:
                indices1[i + 1] = idx1 + 1
            if indices2[i] <= idx1:
                indices2[i] = idx1 + 1
#        print indices1, indices2
#        print indices1-indices2
#        print times

        self.lapIndices = numpy.int_(numpy.vstack([indices1, indices2]).T)

        return self.lapIndices


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
            self.oriented = True
        else:
            print "NOTE: Track already oriented."


    def plot(self, times=[], fig=None, ax=None, offset=2, language='e', chic=False, monochrome=False):

        if not fig:
            fig = pl.figure(figsize=(8, 7))
        if not ax:
            ax = fig.add_subplot(111)

        if not len(times):
            fig, ax, axcb = trajectory.plot(self, fig, ax, offset, language)
            return fig, ax, axcb
        else:
            try:
                self.lapIndices
            except AttributeError:
                sys.exit('No laps calculated! Run traj.getLaps() first!')


        lines = []
        thresh_dx = .5  # in the same units as self.spaceUnit
        for index_pair in self.lapIndices:
            line = []
            for k in numpy.arange(index_pair[0], index_pair[1]):
                if numpy.linalg.norm(self.places[k + 1, 0:2] - self.places[k, 0:2]) > thresh_dx:
                    line.append(tuple(self.places[k, 0:2].tolist()))
                    line.append(tuple((numpy.ones(2) * numpy.nan).tolist()))
                    line.append(tuple(self.places[k + 1, 0:2].tolist()))
                else:
                    line.append(tuple(self.places[k, 0:2].tolist()))
                    line.append(tuple(self.places[k + 1, 0:2].tolist()))
            lines.append(tuple(line))

        # print len(lines), numpy.array(lines).shape, numpy.array(lines)[0]

        if monochrome:
            line_segments = LineCollection(lines, linestyles='solid', colors=numpy.ones(3) * .4, linewidth=1)
        else:
            line_segments = LineCollection(lines, linestyles='solid')
            werte = numpy.arange(times.size + 2)
            line_segments.set_array(werte)
        ax.add_collection(line_segments)

        if not monochrome:
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
            if not monochrome:
                axcb.set_label('Lauf #')
            ax.set_xlabel('x-Position (' + self.spaceUnit + ')')
            ax.set_ylabel('y-Position (' + self.spaceUnit + ')')
        else:  # in english
            if not monochrome:
                axcb.set_label('Run #')
            ax.set_xlabel('x position (' + self.spaceUnit + ')')
            ax.set_ylabel('y position (' + self.spaceUnit + ')')

        # customize plot further?
        if chic:
            custom_plot.allOff(ax)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')

        if not fig:
            pl.show()

        if not monochrome:
            return fig, ax, axcb
        else:
            return fig, ax

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

    def __init__(self, traj, meta={}):
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

    def __init__(self, traj, meta={}):
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

    def __init__(self, traj, meta={}):
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

    def __init__(self, traj, meta={}):
        rewardTrajectory.__init__(self, traj, meta)
        VRlinearMazeTrajectory.__init__(self, self.places, meta, initialized=True)
        # self.times=traj[:, 0]


class decisionMazeRewardTrajectory(rewardTrajectory):

    def __init__(self, traj, meta={}):
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


    def getIRI(self, traj):  # get Interreward interval
        """Calculate the inter-reward interval.

        Due to the teleportation, the inter-reward interval has to be determined
        from the initial position to the end of the lap, i.e., reward time.

        Parameters:
        traj ... the corresponding trajectory
        """

        try:
            traj.lapIndices
        except AttributeError:
            sys.exit('No laps calculated! Run traj.getLaps() first!')

        ici = numpy.array([])
        for index_pair in traj.lapIndices:
            ici = numpy.append(ici, numpy.diff(traj.times[index_pair]))

        return ici


    def getIRD(self, traj):  # get Interreward distance
        """Get the inter-reward distance

        Parameters:
        traj ... the corresponding trajectory

        Returns:
        distances, cumDistances

        """

        distances = numpy.array([])
        cumDistances = numpy.array([])
        for index_pair in traj.lapIndices:
            d = traj.cumPlaces(index_pair[0], index_pair[1])
            cumDistances = numpy.append(cumDistances, d)
            if d.size:
                distances = numpy.append(distances, d[-1])
        return distances, cumDistances


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

