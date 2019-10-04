'''
Created on Jul 9, 2014

@author: chenani
'''
import os, sys
def KKrunner(path,files):
    ''' A function to run KlustaKwik
    This function runs Klustakwik spike sorting program with given parameters!
    
    Parameters
    ----------
    path 
        
    files
        
    '''
    for item in files:
        
        filebase = item.split('.')[0] + ' '
        shank_num = '1 '
        #Number of features (including time)
        # Replace './KlustaKwik' with the path on your system pointing to the executible KlustaKwik
        # DropLastNFeatures 1 means all features except the last will be used
        os.chdir(path)
        os.system(
                  'KlustaKwik '
                  + filebase  
                  + shank_num 
                  +
                  '-UseFeatures' 
                  '-DropLastNFeatures   0'
                  '-UseDistributional   1'
                  '-MaskStarts  0'
                  '-DropLastNFeatures 0 '
                  #'-MaskStarts 300 '
                  '-MaxPossibleClusters 30 '
                  '-PenaltyK   1.0 '
                  '-PenaltyKLogN   0.0 '
                  '-SplitFirst   20 '
                  '-SplitEvery   50'
                  #'-UseDistributional 1 '
                  '-nStarts 3' #The algMaxClustersorithm will be started n times for each inital cluster count between MinClusters and MaxClusters!
                  '-MinClusters   2'
                  '-MaxClusters   20'
                  '-MaxIter 1000'
                  )
def KKrunner2(path,tetrodes):
    for item in tetrodes:
        
        filebase = item + ' '
        shank_num = '1 '
        #Number of features (including time)
        # Replace './KlustaKwik' with the path on your system pointing to the executible KlustaKwik
        # DropLastNFeatures 1 means all features except the last will be used
        os.chdir(path)
        os.system(
                  'KlustaKwik '
                  + filebase  
                  + shank_num 
                  +
                  #'-UseFeatures' 
                  #'-DropLastNFeatures   0'
                  #'-UseDistributional   0'
                  #'-MaskStarts  0'
                  #'-DropLastNFeatures 0 '
                  '-MaxIter 1000'
                  #'-MaskStarts 300 '
                  '-MaxPossibleClusters 30 '
                  #'-PenaltyK   1.0 '
                  #'-PenaltyKLogN   0.0 '
                  '-SplitFirst   20 '
                  '-SplitEvery   50'
                  #'-UseDistributional 1 '
                  #'-nStarts 3' #The algMaxClustersorithm will be started n times for each inital cluster count between MinClusters and MaxClusters!
                  '-MinClusters   2'
                  '-MaxClusters   20'
                  )
