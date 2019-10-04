'''
A module to sequence stuff!
Created on Sep 15, 2014

@author: chenani
'''
import numpy as np

def subsequence(sub,ref):
    '''
    A function to determine the number of specific sequence repeated within a larger sequence.
    
    Parameters
    ----------
    sub: list or an numpy array
        ddd
    
    ref: list or an numpy array
        fff
    
    
    Returns
    ----------
    idx : int array
     This an array containing the index of elements of sub in the ref! if idx is stricktly increasing the sequence sub
    is replayed whithin sequence ref!
    '''
    ref = np.array(ref)
    idx = np.array([])
    for ii in range(len(sub)):
        index = np.where(ref == sub[ii])[0]
        if index.size:
            idx = np.append(idx,index)
        else: 
            return [] , False
    idx = np.int0(idx)
    return idx, True

def median(arr):
    '''mY median... ;)
    A function to pick the median element of an array!
    
    Parameters
    ----------
    arr : numpy.int8
    
    Returns
    ----------
    median : int8
        element of the array corresponding to median!
    '''
    arr = np.array(arr,dtype=np.int8)
    if np.mod(arr.size,2) == 0 :
        return arr[arr.size / 2 - 1]
    else: 
        return arr[arr.size / 2 ]

def sequencer(arr):
    '''sequences the given array!
    A function to sequence a given array with respect to a sequencing rule. i.e. first spike, median spike, center of mass etc.
    
    Parameters
    ----------
    arr : numpy.int8
    '''
    arr = np.array(arr,dtype=np.int8)
    idx_seq = np.array([])
    for item in set(arr):
        idx_seq = np.append(idx_seq,median(np.where(arr == item)[0]))
    idx_seq.sort()

    return arr[np.int0(idx_seq)]

def duplicate_indicator(arr,weights=np.array([])):
    '''Finds and counts the duplicates!
    
    A function to find and count the duplicates in a array of arrays!
    
    '''
    
    arr_cp = np.copy(arr)   
    for ii in range(len(arr)):
        for jj in range(len(arr)):
            if np.array_equal(arr_cp[jj], arr_cp[ii]) and ii != jj:
                arr_cp[jj] = np.array([-1])
    arr_set = []
    for item in arr_cp:
        if item.sum() > 0:
            arr_set.append(item)
    arr_set = np.array(arr_set)
    ### Counting number of repititions
    arr_weights = np.zeros(len(arr_set))
    for ii in range(len(arr_set)):
        for jj in range(len(arr)):
            if np.array_equal(arr[jj],arr_set[ii]):
                if weights.size:
                    arr_weights[ii] += weights[jj]
                else:
                    arr_weights[ii] += 1
    return arr_set,arr_weights
