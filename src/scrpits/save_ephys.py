

'''
A module for reading the Neuralynx files and save as ephys.
written by A. Chenani July 2014
'''
from __future__ import division
import os,sys
sys.path.append("/mnt/Data/Workspaces/Eclipse/dataAnalysis/Sleep-current/src/")
import numpy as np
#import matplotlib as mpl
#mpl.matplotlib_fname()
#import matplotlib.pyplot as pl
import Recordings
from sws_rem_compare import *
import fnmatch



animalPath = "/home/chenani/dataWork/ali/Gerbils/gerbil5/"




tree = []
for (path, dirs, files) in os.walk(animalPath):
    if fnmatch.fnmatch(path,'*sleep*'):
        tree.append(path)     


#Additional Modules
import signale, trajectory,Recordings
# getting the path of the ncs file.
if len(sys.argv) > 1:
    try: 
        folderName = sys.argv[1]
    except:
        print "Usage:",sys.argv[0], "foo.ncs"; sys.exit(1)
else:
    folderName = raw_input('Enter command line arguments: ').split()[0]
src = os.getcwd()
############################PARAMETERS!!!
#NCS loader
onlyWithTTData=False
useRecommended = False
excludeCSCs = []
###T loader
expType = 'real'
prefix = ''
suffix = ''
noSpeck = False
onlyRunning = False
showHeadDir = False
saveFigs = True
useRecommended = False
TTName = '.t'
# other arguments
for argv in sys.argv[2:]:
    if argv.startswith('prefix:'):
        prefix=argv.split(':')[-1]  # name prefix
    if argv.startswith('suffix:'):
        suffix=argv.split(':')[-1]  # name suffix
    if argv.startswith('TT:'):                      # Tetrode files to load
        TTName = argv.split('TT:')[1].strip('[').strip(']')
        TTName = [s for s in TTName.split(',')]
    if argv == 'noSpeck':
        noSpeck = True
    if argv=='onlyRunning':
        onlyRunning = True    
    if argv=='showHeadDir':
        showHeadDir = True            # show head direction
    if argv=='saveFigs':
        saveFigs = True                # save pics
    if argv=='useRecommended':
        useRecommended = True       # use recommendations from metadata.dat
    if argv.startswith('expType:'):
        expType = argv.split(':')[-1]

# initialize in order to make them globally available
###NCS lolader
#cscID = -1
#cscList = signale.NeuralynxCSCList()
###T loader
#loadedSomething = False
#spikes=[]
#ID=-1
#stList=signale.spikezugList(t_start=None, t_stop=None, dims=[2])
#traj = None
#eventData = None

cwd = os.getcwd()


##########################################FUNCTIONS!!!
def getCscData(folderName):
    '''
    a function to load the ncs files, written by 'KT'
    dir sorting added by achenani, oct 2013
    
    '''
    global cscID,cscList,loadedSomething

    if os.path.isdir(folderName):
        dirList=os.listdir(folderName)
        dirList.sort()
        os.chdir(folderName)
    else:
        dirList = [folderName]
        dirList.sort()

    if onlyWithTTData and not any([item.endswith('.t') for item in dirList]):
        os.chdir(cwd)
        sys.exit('The folders do not contain tetrode data (t files)! Therefore skipping folder!')
    for item in dirList:
        if os.path.isfile(item):
            if item.endswith('.ncs') and not any([item.find(str(s))+1 for s in excludeCSCs]):# or item.endswith('2.ncs'):
                print 'loading Neuralynx data', item , 'from folder: '+folderName
                loadedSomething = True
                csc = signale.load_ncsFile(item, showHeader=False)
                cscID += 1
                cscList.append(cscID, csc)
                cscList.addTags(cscID, file=item, dir=folderName)
            elif item.endswith('.raw'):# or item.endswith('2.ncs'):
                print 'loading RAW data', item , 'from folder: '+folderName
                loadedSomething = True
                cscList = []
                cscList =  signale.load_rawFile(item, exclude=excludeCSCs, showHeader=False)
    
        #elif os.path.isdir(item):
        #    getData(item)
    os.chdir('..')
    
def getTData(folderName):
    '''
    a function to load the t files, written by 'KT'
    
    '''
    global spikes, ID, traj, eventData

    if os.path.isdir(folderName):
        dirList=os.listdir(folderName)
        os.chdir(folderName)
    else:
        dirList = [folderName]
        
    dirList.sort()
    for item in dirList:
        print item
        if os.path.isfile(item):
            if (TTName.__class__ == list and item in TTName) or\
                    (TTName.__class__ == str and item.endswith(suffix+'.t') and item.startswith(prefix)):
                print 'loading', item , 'from folder: '+folderName
                spikes = signale.load_tFile(item, showHeader=False)
                ID += 1
                stList.__setitem__(ID, spikes)
                stList.addTags(ID, file=item, dir=folderName)
            # real
            elif expType == 'real':
                if item.endswith('.nvt'):   ## or item.endswith('2.ncs'):
                    print 'loading', item , 'from folder: '+folderName
                    loadedSomething = True
#                     traj = trajectory.load_nvtFile(item, 'linearMaze', showHeader=False)
#                     HDtraj = traj[1]        # head direction
#                     traj = traj[0]          # trajectory
            # vr
            elif expType == 'vr':
                if item.endswith('.nev'):
                    print 'loading', item , 'from folder: '+folderName
                    eventData = signale.load_nevFile(item, showHeader=False)
                elif item.endswith('.traj') and item.find('position')+1\
                 and not item.find('collisions_position')+1 and not item.find('rewardsVisited_position')+1:
                    print 'loading', item , 'from folder: '+folderName
                    if item.startswith('linearMaze'):
                        traj = trajectory.load_trajectory(item, showHeader=False)
                    else:
                        traj = trajectory.load_trajectory(item, showHeader=False)
        elif os.path.isdir(item):
            getTData(item)
    os.chdir('..')

#########################################tasks
for jtem in tree:
    print jtem
    print '***********************'
    fileName = jtem + '/' + jtem.split('/')[-1]
    #############################
    cscID = -1
    cscList = signale.NeuralynxCSCList()
    loadedSomething = False
    ###################
    spikes=[]
    ID=-1
    stList=signale.spikezugList(t_start=None, t_stop=None, dims=[2])
    traj = None
    eventData = None
    #############################
    if (not os.path.exists(fileName+'.ephys')) and jtem.find('leep') > 0:
        
        os.system('ls -l --sort=extension ' + jtem + '/')
        getCscData(jtem)
        getTData(jtem)
        os.chdir(src)
        if len(cscList)>0:
            stList.timeAxis = cscList[0].timeAxis.copy()
            rec = Recordings.ephys(stList,cscList)
            rec.save(fileName)
            del rec
    else:
        print 'There is an ephys file already there or there is nothing to load!!!'

