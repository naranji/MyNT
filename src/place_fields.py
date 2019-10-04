'''
Created on Dec 18, 2013

@author: chenani
'''
from __future__ import division
import os, sys
# add additional custom paths
extraPaths = ["/home/chenani/pypacks/lib/python2.7/site-packages", \
    "/home/thurley/data/", \
    "/home/haas/packages/lib/python2.6/site-packages",
    "/home/chenani/ATLAS-clone/workspace/Sleep/src/", \
    os.path.join(os.path.abspath(os.path.dirname(__file__)), '../scripts')]
for p in extraPaths:
    if not sys.path.count(p):
        sys.path.insert(1, p)
# Additional Modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import signale, trajectory, Recordings
# getting the path of the working directory.
if len(sys.argv) > 1:
    try: 
        folderName = sys.argv[1]
    except:
        print "Usage:", sys.argv[0], "foo.ncs"; sys.exit(1)
else:
    folderName = raw_input('Enter command line arguments: ').split()[0]
src = os.getcwd()
############################PARAMETERS!!!
expType = 'real'
prefix = ''
suffix = ''
noSpeck = False
onlyRunning = False
showHeadDir = False
saveFigs = True
useRecommended = False
TTName = '.t'
###Flags!!!
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
# Initializing spikeList
spikes = []
ID = -1
stList = signale.spikezugList(t_start=None, t_stop=None, dims=[2])
traj = None
eventData = None

cwd = os.getcwd()

#############FUNCTIONS!
def getTData(folderName):
    '''
    a function to load the t files, written by 'KT'
    
    '''
    global spikes, ID, traj, eventData

    if os.path.isdir(folderName):
        dirList = os.listdir(folderName)
        os.chdir(folderName)
    else:
        dirList = [folderName]
        
    dirList.sort()
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
#                    HDtraj = traj[1]        # head direction
                    traj = traj[0]  # trajectory
#         elif os.path.isdir(item):
#             getTData(item)
    os.chdir('..')
    return dirList
###################################################### load data


if os.path.isdir(folderName):
    getTData(folderName)
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

################Setting up trajectories
time, speed = traj.getSpeed()
traj.threshspeed = speed.mean() + 2 * speed.std()  # minimum speed for trajectory as running (m/s)

traj.getTrajDimensions()
traj.getRunningTraj()  # remove from traj.times and
                         # traj.places, times at which animal didn't run <threshspeed
                         
for zaehler, st in enumerate(stList):
    st.traj = traj

    st.getSpikePlaces()
    if showHeadDir:
        st.getSpikeHeadDirection()
    st.getRunningSpikes()  # remove spikes that occurred when the animal was at rest

#####TASKS!
ll = len(stList)
f, farr = pl.subplots(ll, sharex=True,sharey=True)
g, garr = pl.subplots(ll, sharex=True,sharey=True)
for i in range(ll):
    # Spike places
    ylabel = stList.tags[i]['file']
    yl = len(ylabel)
    stList[i].plotSpikesvsPlace(onlyRunning=False, fig=f, ax=farr[i])
    farr[i].set_xlabel('')
    farr[i].set_ylabel (ylabel[yl - 4:yl - 2], fontsize=8)
    farr[i].axes.yaxis.set_ticklabels([])
    lpath = len(stList.tags[0]['dir'])
    f.suptitle(stList.tags[0]['dir'][lpath - 15:lpath] + '\n' + stList.tags[0]['file'][0:6])
    # Field plots
    stList[i].plotField(0.0125, onlyRunning=True, fig=g, ax=garr[i])
    garr[i].set_xlabel('')
    garr[i].set_ylabel (ylabel[yl - 4:yl - 2], fontsize=8)
    garr[i].axes.yaxis.set_ticklabels([])
    g.suptitle(stList.tags[0]['dir'][lpath - 15:lpath] + '\n' + stList.tags[0]['file'][0:6])
    
if saveFigs:
        f.savefig(stList.tags[0]['dir'] + prefix + '_places.png', format='png')
        g.savefig(stList.tags[0]['dir'] + prefix + '_fields.png', format='png')
pl.show()
os.chdir(cwd)
