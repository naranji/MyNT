"""
Show neuralynx ncs-files with python/numpy.
"""
__author__ = ("KT")
__version__ = "2.1, July 2013"



import sys, os, inspect, struct

# add additional custom paths
extraPaths = ["/home/thurley/python/lib/python2.5/site-packages/", \
              "/home/thurley/python/lib/python2.6/site-packages/", \
              "/home/thurley/python/lib/python2.7/dist-packages/", \
              "/home/thurley/data/", \
    os.path.join(os.path.abspath(os.path.dirname(__file__)), '../scripts')]
for p in extraPaths:
    if not sys.path.count(p):
        sys.path.insert(1, p)


# other modules
import numpy

import custom_plot, signale

###################################################### plotting initialization
import matplotlib
matplotlib.matplotlib_fname()
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
import matplotlib as mpl

colors = custom_plot.colors
grau = numpy.array([1, 1, 1]) * .6
grau2 = numpy.array([1, 1, 1]) * .95

Bildformat = 'pdf'

fontsize = 20.0



###################################################### initialization




dummy = sys.argv[1]  # second argument should be the name of the file to load
folderName = dummy.split('\\')[0]
for d in dummy.split('\\')[1:]:
    folderName += '/' + d





# parameters
lang = 'e'
color1 = grau
color2 = 'k'
showFigs = True
saveFigs = False
saveAna = True
onlyWithTTData = False
useRecommended = False
excludeCSCs = []
for argv in sys.argv[2:]:
    if argv == 'noShow':
        showFigs = False  # show pics
    if argv == 'saveFigs':
        saveFigs = True  # save pics
    if argv == 'e' or argv == 'd':  # get language
        lang = argv
    if argv == 'colored':  # colored plots?
        color1 = 'b'
        color2 = 'g'
    if argv == 'onlyWithTTData':  # colored plots?
        onlyWithTTData = True
    if argv == 'useRecommended':
        useRecommended = True  # use recommendations from metadata.dat
    if argv.startswith('exclude:'):
        excludeCSCs = argv.split('exclude:')[1]
        exec 'excludeCSCs =' + excludeCSCs

#################Make a help flag for the script!!!
#
####################################################

# initialize in order to make them globally available
cscID = -1
cscList = signale.NeuralynxCSCList()

loadedSomething = False
cwd = os.getcwd()



if useRecommended:
    if os.path.isfile(folderName + 'metadata.dat'):
        print 'Loading metadata:'
        metadata = signale.read_metadata(folderName + 'metadata.dat', showHeader=True)
        if metadata.has_key('excludeCSCs'):
            exec 'excludeCSCs =' + metadata['excludeCSCs']
            print 'Excluding csc data listed in metadata.dat! CSCs:', excludeCSCs
        print
    else:
        print 'NOTE: There is no metadata.dat. Proceeding without instead.'


###################################################### functions

def getData(folderName):
    global cscID, cscList, loadedSomething

    if os.path.isdir(folderName):
        dirList = os.listdir(folderName)
        os.chdir(folderName)
    else:
        dirList = [folderName]

    if onlyWithTTData and not any([item.endswith('.t') for item in dirList]):
        os.chdir(cwd)
        sys.exit('The folders do not contain tetrode data (t files)! Therefore skipping folder!')
    for item in dirList:
        if os.path.isfile(item):
            if item.endswith('.ncs') and not any([item.find(str(s)) + 1 for s in excludeCSCs]):  # or item.endswith('2.ncs'):
                print 'loading Neuralynx data', item , 'from folder: ' + folderName
                loadedSomething = True
                csc = signale.load_ncsFile(item, showHeader=True)
                cscID += 1
                cscList.append(cscID, csc)
                cscList.addTags(cscID, file=item, dir=folderName)
            elif item.endswith('.raw'):  # or item.endswith('2.ncs'):
                print 'loading RAW data', item , 'from folder: ' + folderName
                loadedSomething = True
                cscList = []
                cscList = signale.load_rawFile(item, exclude=excludeCSCs, showHeader=False)
        # elif os.path.isdir(item):
        #    getData(item)
    os.chdir('..')


def plotData(csc):
    csc.fft()

    fig = pl.figure(figsize=(12, 7))
    fig = pl.figure(figsize=(12, 7))
    pos1 = [.1, .7, .8, .22]
    ax1 = fig.add_axes(pos1)
    pos2 = pos1
    pos2[1] = .39
    ax2 = fig.add_axes(pos2)
    pos3 = pos1
    pos3[1] = .1
    ax3 = fig.add_axes(pos3)

    ax1.set_title(csc.tags['file'])

    csc.plot(fig, ax1)
    custom_plot.huebschMachen(ax1)

    csc.fft_plot(fig, ax2)
    custom_plot.huebschMachen(ax2)
    ax2.set_xlabel('')
    ax2.set_xticklabels([])

    csc.fft_plot(fig, ax3)
    ax3.set_xlim(0, 15.)
    custom_plot.huebschMachen(ax3)






###################################################### main

if os.path.isfile(folderName):
    # get path name of the file to load
    index = folderName.find(folderName.split('/')[-1])
    path = ''
    if index:
        path = folderName[0:index - 1] + '/'

    # load csc data
    csc = signale.load_ncsFile(folderName, showHeader=True)
    loadedSomething = True

    plotData(csc)
elif os.path.isdir(folderName):
    getData(folderName)
    #cscList.changeTimeUnit('s')
    cscList.removeMean()
    fig, ax = cscList.plot()
    fig_fft, ax = cscList.fft_plot(0, 20)
else:
    sys.exit('Folder or data name does not exist.')





###################################################### finishing

os.chdir(cwd)

if saveFigs:
    fig.savefig(folderName + 'csc.png', format='png')
    fig_fft.savefig(folderName + 'csc_fft.png', format='png')
    print 'saved figures ...'


if not loadedSomething:
    sys.exit('The folders do not contain csc data!')

if showFigs:
     pl.show()
raw_input('Press any key to exit!!!')