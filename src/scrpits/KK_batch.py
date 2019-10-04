'''
Created on Jul 30, 2014

@author: chenani
'''
import sys
sys.path.append('/home/chenani/ATLAS-clone/workspace/Sleep/src/')

from extras.trees import locate
from  extras.KK_runner import KKrunner

animal_path = '/home/chenani/DATA-clone/MECLesion_SleepData/Rat434Lesion'
fet_files =[]

for item in locate('*.fet.1',animal_path):
    path = item[0]
    fet_files.append(item[1])
KKrunner(path,fet_files) 
