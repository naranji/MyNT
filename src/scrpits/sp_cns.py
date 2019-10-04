# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 21:30:40 2013

@author: chenani
"""
from __future__ import division
import os,sys
import numpy as np
import matplotlib as mpl
mpl.matplotlib_fname()
import matplotlib.pyplot as pl
# add additional custom paths
extraPaths=["/home/chenani/pypacks/lib/python2.7/site-packages",\
    "/home/thurley/data/",\
    "/home/haas/packages/lib/python2.6/site-packages",\
    os.path.join(os.path.abspath(os.path.dirname(__file__)),'../scripts')]
for p in extraPaths:
    if not sys.path.count(p):
        sys.path.insert(1, p)
#Additional Modules
import signale, trajectory
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

# initialize in order to make them globally available
#NCS lolader
cscID = -1
cscList = signale.NeuralynxCSCList()
loadedSomething = False
#T loader
spikes=[]
ID=-1
stList=signale.spikezugList(t_start=None, t_stop=None, dims=[2])
traj = None
eventData = None

cwd = os.getcwd()

############################FUNCTIONS!!!
def getCscData(folderName):
    '''
    a function to load the ncs files, written by 'KT'
    dir sorting added by achenani, oct 2013
    
    '''
    global cscID, cscList, loadedSomething

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
                csc = signale.load_ncsFile(item, showHeader=True)
                cscID += 1
                cscList.append(cscID, csc)
                cscList.addTags(cscID, file=item, dir=folderName)
            elif item.endswith('.raw'):# or item.endswith('2.ncs'):
                print 'loading RAW data', item , 'from folder: '+folderName
                loadedSomething = True
                cscList = []
                cscList = signale.load_rawFile(item, exclude=excludeCSCs, showHeader=False)
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
                    traj = trajectory.load_nvtFile(item, 'linearMaze', showHeader=False)
                    HDtraj = traj[1]        # head direction
                    traj = traj[0]          # trajectory
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
#######################PLOTTING!

 
 
     
# #########################TASKS!!!
os.system('ls -l --sort=extension ' + folderName)
getCscData(folderName)
getTData(folderName)
os.chdir(src)
################################################################################
################################################################################
# no = 0
#if len(cscList) % 2:
#    f,axarr = pl.subplots(int(len(cscList)/2) + 1,2)
#else:
#    f,axarr = pl.subplots(int(len(cscList)/2),2)
# for item in cscList:
#     item.filter(150,250)
#     item.hilbertTransform()
#     item.rms(10)
#     item.ripple_recorder()


#    axarr[int(no / 2),no % 2].hist((item.ripples[:,1] -item.ripples[:,0]) * item.dt , bins = 75)
#    axarr[int(no / 2),no % 2].set_title('Channel: '+ str(item.tags.get('channel')), fontsize=10, fontweight='bold')    
#    no+=1
#f.suptitle('Width Distribution in ' + folderName.split('/')[-2])
################################################################################
################################################################################
# stl = stList.time_slice(cscList.t_start,int(cscList.t_start + ((cscList.t_stop-cscList.t_start)/10)))
# stList.timeAxis = cscList[0].timeAxis.copy()
# stList.burstDetector()
# n_pop = np.zeros(30)
# pop_size = np.zeros(n_pop.size)
# for k in range(n_pop.size):
#     stList.burstDetector('gauss', (k+1)*10)
#     n_pop[k] = stList.bursts.shape[0]
#     a = []
#     for r in range(stList.bursts.shape[0]):
#         a.append(stList.bursts[r][1] - stList.bursts[r][0])
#     a = np.array(a)
#     pop_size[k] = a.mean()
#     print '%d pops using %d window'%(n_pop[k],k)

################################################################################
################################################################################
# 
# n_ripp = np.zeros(30)
# ripp_size = np.zeros(n_ripp.size)
# for k in range(n_ripp.size):
#     cscList.rippleDetect((k+1)*10)
#     n_ripp[k] = cscList.ripples.shape[0]
#     a = []
#     for r in range(cscList.ripples.shape[0]):
#         a.append(cscList.ripples[r][1] - cscList.ripples[r][0])
#     a = np.array(a)
#     ripp_size[k] = a.mean()
#     print '%d ripples using %d window'%(n_ripp[k],k)
# f,plarr = pl.subplots(3, 1, 1)
# plarr[0].plot(10 * np.arange(1,n_ripp.size + 1),n_ripp)
# plarr[2].set_xlabel("Sigma")
 



     
