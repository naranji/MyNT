'''
Created on Jul 30, 2014
This is a module to tun KlustaKwik on parallel
@author: chenani
'''
import sys
sys.path.append('/home/chenani/ATLAS-clone/workspace/Sleep/src/'),
sys.path.append( '/home/chenani/Downloads/pprocess-0.5/')
from extras.trees import locate
from  extras.KK_runner import KKrunner2
import pprocess
animal_path = '/home/chenani/DATA-clone/MECLesion_SleepData/Rat434Lesion'
fet_files =[]
path = []
for item in locate('*.fet.1',animal_path):
    path = item[0]
    fet_files.append(item[1])
#KKrunner(path,fet_files) 
#prepare parallelization
list_of_TTs = ['TT1','TT2']
nproc = 4      # maximum number of simultaneous processes desired
results = pprocess.Map(limit=nproc, reuse=1)
parallel_function = results.manage(pprocess.MakeReusable(KKrunner2))
[parallel_function(path,args) for args in list_of_TTs];  # Start computing things
print fet_files