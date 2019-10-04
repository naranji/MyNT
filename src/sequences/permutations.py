'''
A module that helps to permutate stuff!
Created on Sep 15, 2014

@author: chenani
'''
from sequences.seqz import *
import numpy as np
import scipy.misc as scmsc

def next_permutation(arr):
    '''
    Computes the next lexicographical permutation of the specified list in place,
    returning whether a next permutation existed. (Returns False when the argument
    is already the last possible permutation.)
    
    Parameters
    ----------
    arr : list
        This is initial list, which would be modified with the next lexicographical permutation(if exists!!!)
    
    Returns
    -------
    bool
        Tells you whether a next permutation exists or what?! ;)
    
    Example
    -------
        arr = [0, 1, 0]
        next_permutation(arr)  (returns True)
        arr has been modified to be [1, 0, 0]
    
    References
    ----------
    Next lexicographical permutation algorithm (Python)
    By Nayuki Minase, 2014. Public domain.
                 http://nayuki.eigenstate.org/page/next-lexicographical-permutation-algorithm
    
    '''
    # Find non-increasing suffix
    i = len(arr) - 1
    while i > 0 and arr[i - 1] >= arr[i]:
        i -= 1
    if i <= 0:
        return False
    
    # Find successor to pivot
    j = len(arr) - 1
    while arr[j] <= arr[i - 1]:
        j -= 1
    arr[i - 1], arr[j] = arr[j], arr[i - 1]
    
    # Reverse suffix
    arr[i : ] = arr[len(arr) - 1 : i - 1 : -1]
    return True

def next_permutation_comp(arr, comp):
    '''
        Computes the next lexicographical permutation of the specified list in place,
     returning whether a next permutation existed. (Returns False when the argument
     is already the last possible permutation.)
 
     comp is a comparison function - comp(x, y) returns a negative number if x is considered to be less than y,
     a positive number if x is considered to be greater than y, or 0 if x is considered to be equal to y.
 
    References
    ----------
    Next lexicographical permutation algorithm (Python)
    By Nayuki Minase, 2014. Public domain.
                 http://nayuki.eigenstate.org/page/next-lexicographical-permutation-algorithm
    '''
    # Find non-increasing suffix
    i = len(arr) - 1
    while i > 0 and comp(arr[i - 1], arr[i]) >= 0:
        i -= 1
    if i <= 0:
        return False
    
    # Find successor to pivot
    j = len(arr) - 1
    while comp(arr[j], arr[i - 1]) <= 0:
        j -= 1
    arr[i - 1], arr[j] = arr[j], arr[i - 1]
    
    # Reverse suffix
    arr[i : ] = arr[len(arr) - 1 : i - 1 : -1]
    return True

def number_of_permutations(arr):
    Cnk = [] #keep the combinations!
    s = 0    # sum of repitiotions of elements in arr!
    arrcp = np.array(arr).copy()
    arrcp.sort()
    for item in set(arrcp):
        reps = np.where(arrcp == item)[0].size
        Cnk.append(scmsc.comb(arrcp.size - s,reps,True))
        s += reps
    return np.prod(np.array(Cnk))

def all_permutations(arr):
    '''
    produces all possible permutations of a given array using lexographical ordering.
    Sequence the arrays and returns the weights of all possible sequences!
    '''
    sqStack = []
    weightsStack = []
    for item in arr:
        item.sort()
        arr_seqz = [sequencer(item)]
        
        #####creating MUA permutation list!
        while next_permutation(item):
            arr_seqz.append(sequencer(item))
        arr_seqz = np.array(arr_seqz) 
        
        #####Sequencing Permutations
        arr_seqz_set = np.array([])
        sample_seq = arr_seqz[1].copy()
        sample_seq.sort()
        arr_seqz_set = np.append(arr_seqz_set,sample_seq)
    
        while next_permutation(sample_seq):
            arr_seqz_set = np.append(arr_seqz_set,sample_seq)
        arr_seqz_set = arr_seqz_set.reshape(np.math.factorial(sample_seq.size),arr_seqz[0].size)
        
        #####Calculating the weights of all possible sequences
        weights = []
        for item in arr_seqz_set:
            selection = arr_seqz[np.where(arr_seqz[:,0]==item[0])[0]]
            for ii in range(1,arr_seqz_set.shape[1]):
                selection = selection[np.where(selection[:,ii]==item[ii])[0]]
            weights.append(selection.size/float(arr_seqz.size))
        weightsStack.append(np.array(weights))
        sqStack.append(arr_seqz_set)
    return sqStack,weightsStack