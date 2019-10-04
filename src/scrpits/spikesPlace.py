"""
Load (AD Redishs's MClust) t-files with python/numpy and plot spikes
against place.

Use expType:real, for real environment data, i.e., with position tracked into Neuralynx NVT files.
Use expType:vr, for virtual reality data, i.e., with position tracked into traj files.

"""
__author__ = ("KT", "Moritz Dittmeyer")
__version__ = "4.2, July 2013"

# python modules
import sys, os, inspect, struct

# Uncomment this line when you want to set the desired arguments (ex. path of t files)
# after running the script!

# sys.argv[1:] = raw_input('Enter command line arguments: ').split()


# add additional custom paths
extraPaths = ["/home/chenani/pypacks/lib/python2.7/site-packages", \
              "/home/thurley/python/lib/python2.5/site-packages/", \
              "/home/thurley/python/lib/python2.6/site-packages/", \
              "/home/thurley/python/lib/python2.7/dist-packages/", \
              "/home/thurley/data/", \
              "/home/haas/packages/lib/python2.6/site-packages", \
    os.path.join(os.path.abspath(os.path.dirname(__file__)), '../scripts')]
for p in extraPaths:
    if not sys.path.count(p):
        sys.path.insert(1, p)


# other modules
import numpy
import NeuroTools.signals as NTsig

# custom made modules
import signale, trajectory, custom_plot



###################################################### plotting initialization

import matplotlib as mpl
import matplotlib.pyplot as pl

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection, CircleCollection
# Changing the __init__ method of AutoLocator 
def AutoLocatorInit(self):
    mpl.ticker.MaxNLocator.__init__(self, nbins=4)
mpl.ticker.AutoLocator.__init__ = AutoLocatorInit


fontsize = 20.0
markersize = 6

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = markersize
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'


colors = ['#FF0000', '#0000FF', '#008000', '#00FFFF', '#FF00FF', '#EE82EE',
        '#808000', '#800080', '#FF6347', '#FFFF00', '#9ACD32', '#4B0082',
        '#FFFACD', '#C0C0C0', '#A0522D', '#FA8072', '#FFEFD5', '#E6E6FA',
        '#F1FAC1', '#C5C5C5', '#A152ED', '#FADD72', '#F0EFD0', '#EEE6FF',
        '#01FAC1', '#F5F5F5', '#A152FF', '#FAFD72', '#F0EFDF', '#EEEFFF',
        '#F1FA99', '#C9C9C9', '#A152DD', '#FA5572', '#FFFFD0', '#EDD6FF']

grey = numpy.ones(3) * .25
transred = '#FFA1A1'



###################################################### commandline paramters

dummy = sys.argv[1]  # second argument should be the name of the folder to load
#This part is for fixing the path in the windows machines!!!
folderName = dummy.split('\\')[0]
for d in dummy.split('\\')[1:]:
    folderName += '/' + d


# parameters
expType = 'real'
prefix = ''
suffix = ''
noSpeck = False
onlyRunning = False
showHeadDir = False
saveFigs = True
useRecommended = False
TTName = '.t'
for argv in sys.argv[2:]:
    if argv.startswith('prefix:'):
        prefix = argv.split(':')[-1]  # name prefix
    if argv.startswith('suffix:'):
        suffix = argv.split(':')[-1]  # name suffix
    if argv.startswith('TT:'):  # Tetrode files to load
        TTName = argv.split('TT:')[1].strip('[').strip(']')
        TTName = [s for s in TTName.split(',')]
    if argv == 'noSpeck':
        noSpeck = True
    if argv == 'onlyRunning':
        onlyRunning = True	
    if argv == 'showHeadDir':
        showHeadDir = True  # show head direction
    if argv == 'saveFigs':
        saveFigs = True  # save pics
    if argv == 'useRecommended':
        useRecommended = True  # use recommendations from metadata.dat
    if argv.startswith('expType:'):
        expType = argv.split(':')[-1]


###################################################### initialization


# initialize in order to make them available globally
spikes = []
ID = -1
stList = signale.spikezugList(t_start=None, t_stop=None, dims=[2])
traj = None
eventData = None

cwd = os.getcwd()


if useRecommended:
    if os.path.isfile(folderName + 'metadata.dat'):
        print 'Loading metadata:'
        metadata = signale._read_metadata(folderName + 'metadata.dat', showHeader=True)
        if metadata.has_key('tt'):
            exec 'TTName =' + metadata['tt']
            print 'Taking tetrode data listed in metadata.dat! TT:', TTName
        print
        print
    else:
        print 'NOTE: There is no metadata.dat. Proceeding without instead.'


###################################################### functions

def getData(folderName):
    global spikes, ID, traj, eventData

    if os.path.isdir(folderName):
        dirList = os.listdir(folderName)
        os.chdir(folderName)
    else:
        dirList = [folderName]

    for item in dirList:
        if os.path.isfile(item):
            if (TTName.__class__ == list and item in TTName) or\
                    (TTName.__class__ == str and item.endswith(suffix + '.t') and item.startswith(prefix)):
                print 'loading', item , 'from folder: ' + folderName
                spikes = signale.load_tFile(item, showHeader=False)
                ID += 1
                stList.__setitem__(ID, spikes)
                stList.addTags(ID, file=item, dir=folderName)
            # real
            elif expType == 'real':
                if item.endswith('.nvt'):  # # or item.endswith('2.ncs'):
                    print 'loading', item , 'from folder: ' + folderName
                    loadedSomething = True
                    traj = trajectory.load_nvtFile(item, 'linearMaze', showHeader=False)
                    HDtraj = traj[1]  # head direction
                    traj = traj[0]  # trajectory
            # vr
            elif expType == 'vr':
                if item.endswith('.nev'):
                    print 'loading', item , 'from folder: ' + folderName
                    eventData = signale.load_nevFile(item, showHeader=False)
                elif item.endswith('.traj') and item.find('position') + 1\
                 and not item.find('collisions_position') + 1 and not item.find('rewardsVisited_position') + 1:
                    print 'loading', item , 'from folder: ' + folderName
                    if item.startswith('linearMaze'):
                        traj = trajectory.load_trajectory(item, showHeader=False)
                    else:
                        traj = trajectory.load_trajectory(item, showHeader=False)
        elif os.path.isdir(item):
            getData(item)
    os.chdir('..')


###################################################### load data


if os.path.isdir(folderName):
    getData(folderName)
elif os.path.isfile(folderName):
    sys.exit('Point to a folder not a single file.')
else:
    sys.exit('Folder or data name does not exist.')
os.chdir(cwd)

if not ID + 1:
    sys.exit('The folders do not contain spike data!')

# real
elif expType == 'real':
    # set begining to zero time
    # for st in stList:
    #    st.spike_times -= st.t_start
    # stList._spikezugList__recalc_startstop()

    # cut away all earlier spikes
    # stList = stList.time_slice(0., stList.t_stop)

    # change to seconds since trajectories are also stored in seconds
    stList._spikezugList__recalc_startstop()
    stList.changeTimeUnit('s')

    # stList = stList.time_slice(0, 20000)      # reduce data size a bit

# vr
elif expType == 'vr':
    # Remove the first two entries, when the data was recorded in the VR setup
    #       since they are just initialization.
    eventData.times = eventData.times[2:]
    eventData.eventStrings = eventData.eventStrings[2:]
    eventData.t_start = eventData.times[0]


    # set begining to zero time
# #    eventData.times -= eventData.t_start
# #    for st in stList:
# #        st.spike_times -= eventData.t_start
# #    stList._spikezugList__recalc_startstop()


    # cut away all earlier spikes
    # stList = stList.time_slice(0., stList.t_stop)

    # change to seconds since trajectories are also stored in seconds
    eventData.changeTimeUnit('s')
    stList.changeTimeUnit('s')

    traj.times = eventData.times[:-2]  # ordentlich machen!
    if traj.places.shape[0] > traj.times.shape[0]:  # nur ganz schlechter workaround, 24.03.13!
        traj.places = traj.places[:traj.times.shape[0], :]
    elif traj.places.shape[0] < traj.times.shape[0]:
        traj.times = traj.times[:traj.places.shape[0]]
# #    traj.times += eventData.t_start

    # stList = stList.time_slice(0, 20000)      # reduce data size a bit



###################################################### main

###################################################### clean data
time, speed = traj.getSpeed()
traj.threshspeed = speed.mean() + 2 * speed.std()  # minimum speed for trajectory as running (m/s)

traj.getTrajDimensions()
traj.getRunningTraj()  # remove from traj.times and
                         # traj.places, times at which animal didn't run <threshspeed


###################################################### spikes against place

for zaehler, st in enumerate(stList):
    st.traj = traj

    st.getSpikePlaces()
    if showHeadDir:
        st.getSpikeHeadDirection()

# #    # shorten traj to test
# #    traj.times=traj.times[0:len(traj.times)/10]
# #    traj.places=traj.places[0:len(traj.places)/10]

    st.getRunningSpikes()  # remove spikes that occurred when the animal was at rest


###################################################### plots



    if noSpeck:
        fig = pl.figure(20 + zaehler, figsize=(12, 6))
        numRows = 2
    else:
        fig = pl.figure(20 + zaehler, figsize=(12, 8))
        numRows = 4


    ######################### spike/place -> running & not running in color transred

    ax = fig.add_subplot(numRows, 1, 1)

    if not noSpeck:
        ax.set_position([ 0.125, 0.79, 0.775, 0.17])

    st.plotSpikesvsPlace(onlyRunning=onlyRunning, showHeadDir=showHeadDir, fig=fig, ax=ax)

    custom_plot.huebschMachen(ax)



    ########################################################## spike/place -> PCOLOR_PLOT

    if not noSpeck:
        ax = fig.add_subplot(numRows, 1, 2)

        st.plotField(.03, onlyRunning=True, fig=fig, ax=ax)

        ax.set_xlabel('', visible=False)
        ax.set_ylabel('', visible=False)

        ax.set_xlim(st.traj.xlim + numpy.diff(st.traj.xlim) * .05 * [-1, 1])
        ax.set_ylim(st.traj.ylim + numpy.diff(st.traj.ylim) * .1 * [-1, 1])
        custom_plot.allOff(ax)


    ##################################    spike/time plots


    ########################################### 1

    if noSpeck:
        ax = fig.add_subplot(numRows, 1, 2)
    else:
        ax = fig.add_subplot(numRows, 1, 3)

    ax.plot(st.traj.run_times, st.traj.run_places[:, 0], linewidth=2.5, color=grey)  # only run trajectory
    ax.plot(st.traj.times, st.traj.places[:, 0], linewidth=.5, color='k')  # original/full trajectory
    ax.plot(st.spike_times, st.spikePlaces[:, 0], color=transred, marker='.', linestyle='None')  # all spiketimes
    ax.plot(st.run_spikeTimes, st.run_spikePlaces[:, 0], 'r.')  # only running-spikes

    if noSpeck:
        ax.set_xlabel('Time (s)')
    else:
        ax.set_xticklabels([])
    ax.set_ylabel('x position (' + st.traj.spaceUnit + ')')
    ax.set_xlim(traj.times[0], st.traj.times[-1])
    custom_plot.huebschMachen(ax)

    ############### running-right spike/time

    if not noSpeck:
        ax = fig.add_subplot(numRows, 1, 4)

        st.traj.getLeftAndRightwardRuns(onlyRunning=onlyRunning)
        st.getLeftAndRightwardSpikes(onlyRunning=onlyRunning)

        # plot all spikes! transred .'s
        ax.plot(st.spike_times, st.spikePlaces[:, 0], color=transred, marker='.', linestyle='None')

        ax.plot(st.traj.rechts_times, st.traj.rechts_places[:, 0], linewidth=1.0, color='m', label="right")  # plot traj running right
        if hasattr(st, 'rechts_spikeTimes'):
            ax.plot(st.rechts_spikeTimes, st.rechts_spikePlaces[:, 0], color='r', marker='o', markersize=4, linestyle='None')  # plot spikes right red o's


    ########running-left spike/time#######

        ax.plot(st.traj.links_times, st.traj.links_places[:, 0], linewidth=1.0, color='b', label="left")  # plot traj running left
        if hasattr(st, 'links_spikeTimes'):
            ax.plot(st.links_spikeTimes, st.links_spikePlaces[:, 0], color='g', marker='o', markersize=4, linestyle='None')  # plot spikes left green o's

        ax.set_xlabel('Time (' + st.traj.timeUnit + ')')
        ax.set_xlim(traj.times[0], st.traj.times[-1])

        custom_plot.huebschMachen(ax)


    #--- text for data set
    # .995, .98
    fig.text(.995, .01, stList.tags[zaehler]['dir'] + ', ' + stList.tags[zaehler]['file']\
        + ', v_run >= ' + str(traj.threshspeed) + ' ' + traj.spaceUnit + '/' + traj.timeUnit, \
        fontsize=fontsize - 10, horizontalalignment='right')

    if noSpeck:
        fig.tight_layout()

    if saveFigs:
        fig.savefig(stList.tags[zaehler]['dir'] + stList.tags[zaehler]['file'].split('.')[0] + '.png', \
            format='png')  # save figure , dpi=300


###################################################### spikes against speed


# #time, speed = traj.getSpeed()
# #
# #numSpikes = []
# #
# #for i in range(time.shape[0]-1):
# #    numSpikes.append(stList.spiketrains[0].time_slice(time[i], time[i+1]).spike_times.shape[0])
# #
# #numSpikes = numpy.array(numSpikes)
# #
# #fig = pl.figure(21)
# #ax = fig.add_subplot(121)
# #ax.plot(time, speed)
# #ax.plot(time[:-1], numSpikes, '.')
# #ax.set_xlabel('Time (s)')
# #ax.set_ylabel('Speed (m/s) / Number of spikes')
# #
# #ax = fig.add_subplot(122)
# #ax.plot(speed[:-1], numSpikes, '.')
# #ax.set_xlabel('Speed (m/s)')
# #ax.set_ylabel('Number of spikes')
# #
# #
# #



###################################################### finishing



pl.show()



