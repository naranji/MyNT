"""
trajectory.io
=============

A module for input/output functions for trajectories.
"""
__author__ = ("KT")
__version__ = "7.0, August 2013"


# python modules
import struct, os

# other modules
import numpy

# custom made modules

# package modules
import vr
import trajectory

###################################################### FUNCTIONS


def load_trajectory(fileName, showHeader=False):
    """
    For loading VR trajectory data
    with python/numpy. Returns a trajectory object.
    """

    traj = numpy.loadtxt(fileName)

    if traj.shape.__len__() == 1:
        traj = traj.reshape(1, traj.shape[-1])
    meta = _read_metadata(fileName, showHeader=showHeader)

    # in case there are no values in the file
    if traj.size < 5 and traj.mean() == 0:
        traj = numpy.array([[0., 0., 0., 0.]])

    if not meta.has_key('mazetype'):
        meta['mazetype'] = ''


    if meta['mazetype'].find('linearMaze') + 1:
        return vr.VRlinearMazeTrajectory(traj, meta)
    elif meta['mazetype'].find('yDecision') + 1:
        return vr.decisionMazeTrajectory(traj, meta)
    else:
        return vr.VRtrajectory(traj, meta)


def load_collisionTrajectory(fileName):

    traj = numpy.loadtxt(fileName)
    if traj.shape.__len__() == 1:
        traj = traj.reshape(1, traj.shape[-1])

    # in case there are no values in the file
    if traj.size < 5 and traj.mean() == 0:
        traj = numpy.array([[0., 0., 0., 0.]])

    meta = _read_metadata(fileName)
    if meta['mazetype'] == 'linearMaze':
        return vr.linearMazeCollisionTrajectory(traj, meta)
    else:
        return vr.collisionTrajectory(traj, meta)



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

        count = -1
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
                if count < 0 and showHeader:
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
    if 'VT_parameters.dat' in os.listdir(filePath):  # load video tracker parameters
        meta.update(_read_metadata(filePath + '/VT_parameters.dat'))

    if meta.has_key('mazetype'):
        mazeType == meta['mazetype']

    if mazeType == 'linearMaze':
        meta.update({'mazetype' : mazeType})
        posTraj = trajectory.linearMazeTrajectory(posSamples, meta)
        eulTraj = []
        # eulTraj = linearMazeTrajectory(eulSamples, meta)
    else:
        posTraj = trajectory.trajectory(posSamples, meta)
        eulTraj = trajectory.trajectory(eulSamples, meta)

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
        return vr.linearMazeRewardTrajectory(traj, meta)
    elif meta['mazetype'] in ['teleporter']:
        return vr.rewardTrajectory(traj, meta)
    elif meta['mazetype'].find('yDecision') + 1:
        return vr.decisionMazeRewardTrajectory(traj, meta)
    else:
        return vr.rewardTrajectory(traj, meta)

# auskommentiert am 07.05.13, da unklar warum man das hier drin braucht
# #def load_rewardsPos(fileName):
# #    """
# #    For loading reward related information for a maze coming from Blender.
# #    """
# #
# #    # initialize
# #    rewardsPosNames = []
# #    rewardsPosNamesAdditionals = {}
# #    rewardsPos = []
# #    initialPos = []
# #    rewardsArea = []
# #    rewardsAreaType = ''
# #
# #    meta = _read_metadata(fileName)
# #    for key in meta:
# #        if key=='rewardsArea':
# #            exec 'dummy='+str(meta[key])
# #            if not isinstance(dummy, list):
# #                dummy = [dummy]
# #            rewardsArea.extend(dummy)
# #        elif key=='rewardsRadius':
# #            rewardsArea=[meta[key]]
# #            print 'NOTE: loading reward positions with deprecated header style using rewardsRadius'
# #        elif key=='rewardsAreaType':
# #            rewardsAreaType=meta['rewardsAreaType']
# #        elif key.endswith('add'):
# #            exec 'dummy='+meta[key]
# #            rewardsPosNamesAdditionals[key.split('add')[0]]=dummy
# #        elif key == 'initialPos':
# #            exec 'dummy='+meta[key]
# #            initialPos.extend(dummy)
# #        else:
# #            rewardsPosNames.append(key)
# #            exec 'dummy='+meta[key]
# #            rewardsPos.append(dummy)
# #    if not len(initialPos):
# #        initialPos = rewardsPos[0]
# #        print 'NOTE: loaded obsolete maze style without initialPos.'
# #
# #    return rewardsPosNames, rewardsPosNamesAdditionals, rewardsPos, rewardsArea, rewardsAreaType, initialPos


def load_params(fileName):
    """ For loading the parameters from tsv parameters files."""

    params = numpy.loadtxt(fileName, delimiter='\t', dtype='str')  # tsv fromat => tab delimiters
    if params.shape.__len__() == 1:
        params = numpy.array([params])

    parameters = []
    times = []
    for p in params:
        dictio = {}
        t = []
        # print p
        for st in p:
            # print st
            if st.split(':').__len__() > 1:
                exec 'dictio[\'' + st.split(':')[0] + '\']=' + st.split(':')[1]
            else:
                t = numpy.append(t, st)  # NOTE: not exported so far
        dictio['scaling'] = numpy.array(dictio['scaling'])
        times.append(t)
        parameters.append(dictio)
    times = numpy.float_(times).flatten()

    meta = _read_metadata(fileName)

    return trajectory.paramsTrajectory(times, parameters, meta)


def load_stimuli(stimuliDates, folderNames, rewardTraj, folderName):
    """For loading the cvs stimuli files for the yDecisionwCues experiments.

        Rather specifically tuned to the yDecisionwCues experiments.
        I don't really like it to be in here, but need to use it several times
        and don't have a better place for now.
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
        # print stimuliFile, stimuliInFile
        print 'file:', d
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
        # print 'file:', stimuliDates[i]
        # print 'presentations of sets 0-3 (%):', numpy.round(numpy.sum(decisions_sum_perSubSet[:-2])/f.shape[0], 2)
        # print 'presentations of set 4 (%):', numpy.round(numpy.sum(decisions_sum_perSubSet[8:])/f.shape[0], 2)
        # print 'number of stimuli in file:', f.shape[0]


    return stimuli, stimulusSetsNames, stimulusSets



def _read_metadata(fileName, showHeader=False):
    """
    Read the informations that may be contained in the header of
    the trajectory object, if saved in a text file.
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

    metadata['file'] = os.path.abspath(fileName)

    return metadata


###################################################### CLASSES
