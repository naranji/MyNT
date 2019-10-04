'''
Created on Jun 25, 2014

@author: chenani
'''
from __future__ import division
import os,sys
import numpy as np
import pickle as pkl
import fnmatch
import os.path
import collections
##########################################
dataFolder = "/home/chenani/DATA-clone/"
animalPath = "/home/chenani/DATA-clone/"
##########################################FUNCTIONS
def path_finder(data=dataFolder):
    '''Locate all paths for different animals in a data folder.
    '''
    paths = []
    for root, dirs,files in os.walk(os.path.abspath(data)):
        dirs[:] = [d for d in dirs if d.find('erbil') > 0]
        dirs = [os.path.join(root,d) for d in dirs]
        [paths.append(d) for d in dirs]
    return paths        
def file_locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
        supplied root directory.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield [path,filename]
def save(report, filename):
        '''
        This function saves!!
        '''
        pkl.dump(report, open(filename + ".rpt", 'wb'), pkl.HIGHEST_PROTOCOL)
        print "data has been saved to %s.rpt" % filename
##########################################
animalList = path_finder(animalPath)
for animal in animalList:
    vr = animal+'/vr'
    real = animal+'/real'
    real_rpts = []
    vr_rpts = []
    for session in file_locate('*.rpt',real):
        real_rpts.append(os.path.join(session[0],session[1]))
    for session in file_locate('*.rpt',vr):
        vr_rpts.append(os.path.join(session[0],session[1]))
        
    report = collections.OrderedDict()
    report['animal_name'] = animal.split('/')[-1]
    s_sum = 0
    a_sum = 0
    sleep_duration = 0
    number_of_ripples = 0
    rem_content = 0
    for item in vr_rpts:
        if item.find('leep') > 1:
            s_sum+=1
            rpt = pkl.load(open(item, "rb"))
            sleep_duration += rpt['duration (min)']
        else:
            a_sum+=1
            
    report['number of VR behaving sessions'] = a_sum
    report['number of VR sleep sessions'] = s_sum
    s_sum = 0
    a_sum = 0
    for item in real_rpts:
        if item.find('leep') > 1:
            s_sum+=1
        else:
            a_sum+=1
    report['number of REAL behaving sessions'] = a_sum
    report['number of REAL sleep sessions'] = s_sum
    print report


#########################################
# for item in rptList:
#     print item
#     if str(item).find('leep') > 1:
#         print os.path.join(item[0],item[1])
#         rptFile = os.path.join(item[0],item[1])
#         report = pkl.load(open(rptFile, "rb"))
#     