'''
A module to copy rpt files.
Created on Aug 9, 2014

@author: chenani
'''
import os
import pickle
from trees import locate

animalPath = '/home/chenani/DATA-clone/Gerbil_9588'
rptList = []
for item in locate("*.rpt",animalPath):
    rptList.append(os.path.join(item[0],item[1]))
######################################################
destFolder = animalPath + '/test'
mkdir_cmd = 'mkdir ' + destFolder
mkdir1_cmd = 'mkdir ' + destFolder + '/vr'
mkdir2_cmd = 'mkdir ' + destFolder + '/real'
os.system(mkdir_cmd)
os.system(mkdir1_cmd)
os.system(mkdir2_cmd)
#######################
i = 1
for item in rptList:
    if item.find('real') > 0:
        if os.path.exists( destFolder + '/real/report.rpt'):
            os.system('mv ' + destFolder + '/real/report.rpt ' + destFolder + '/real/report'+ str(i) + '.rpt' )
            i+=1
        os.system('cp ' + item + ' '  + animalPath + '/test/real')
    elif item.find('vr') > 0:
        if os.path.exists( destFolder + '/vr/report.rpt'):
            os.system('mv ' + destFolder + '/vr/report.rpt ' + destFolder + '/vr/report'+ str(i) + '.rpt' )
            i+=1
        os.system('cp ' + item + ' '  + animalPath + '/test/vr')
i = 0        