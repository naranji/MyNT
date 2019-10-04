from scipy.misc import comb
import scipy as scp
import numpy as np
import matplotlib.pyplot as pl
import pickle as pkl
import sets
import random
import itertools
import os,sys,fnmatch
import timeit
import scipy.stats


# FUNCTIONS

def subsequence(sub,ref):
    '''
    A function to determine the number of specefic sequence repeated within a larger sequence.
    
    Parameters:
    sub:
    ref:
    
    
    Returns:
    idx : This an array containing the index of elements of sub in the ref! if idx is stricktly increasing the sequence sub
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
    '''
    mY median... ;)
    '''
    arr = np.array(arr)
    if np.mod(arr.size,2) == 0 :
        return arr[arr.size / 2 - 1]
    else: 
        return arr[arr.size / 2 ]
    
def sequencer(arr,method = 'median'):
    '''
    This function sequences the given array(with possible repeated elements) into an array of distinct elements! Considering either 
    first or the midian position of repeated elements.
    '''
    arr = np.array(arr)
    idx_seq = np.array([])
    if method == 'median':
        for item in set(arr):
            idx_seq = np.append(idx_seq,median(np.where(arr == item)[0]))
        idx_seq.sort()
    if method == 'first':
        for item in set(arr):
            idx_seq = np.append(idx_seq,np.where(arr == item)[0][0])
        idx_seq.sort()

    return arr[np.int0(idx_seq)]

def duplicate_indicator(arr_group,weights=np.array([])):
    '''
    Finds and counts the duplicates! This Function makes a frequency distribution the number of
    arrays within a goup of arrays.
    
    Parameters:
    -----------
    arr_group: The one that you want to count
    weights: an array containig the wights(results of a previous counting of etc.). 
    If non empty the distribution will be weighted 
             with respect to this array!
    Returns:
    -----------
    arr_set: Set of distinct arrays in arr_group
    arr_weights: repitition counts of elements in arr_set.
    '''
    
    arr_cp = np.copy(arr_group)   
    for ii in range(len(arr_group)):
        for jj in range(len(arr_group)):
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
        for jj in range(len(arr_group)):
            if np.array_equal(arr_group[jj],arr_set[ii]):
                if weights.size:
                    arr_weights[ii] += weights[jj]
                else:
                    arr_weights[ii] += 1
    return arr_set,arr_weights

def next_permutation(arr):
    '''
    
        Computes the next lexicographical permutation of the specified list in place,
     returning whether a next permutation existed. (Returns False when the argument
     is already the last possible permutation.)
    
    
         Example:
            arr = [0, 1, 0]
            next_permutation(arr)  (returns True)
            arr has been modified to be [1, 0, 0]
    Reference:
    -----------
    Nayuki Minase, 2014. Public domain.
    http://nayuki.eigenstate.org/page/next-lexicographical-permutation-algorithm
            '''
    
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

def number_of_permutations(arr):
    Cnk = [] #keep the combinations!
    s = 0    # sum of repitiotions of elements in arr!
    arrcp = np.array(arr).copy()
    arrcp.sort()
    for item in set(arrcp):
        reps = np.where(arrcp == item)[0].size
        Cnk.append(comb(arrcp.size - s,reps,True))
        s += reps
    return np.prod(np.array(Cnk))
def all_permutations(arr):
    '''
    produces all possible permutations of a given array using lexographical ordering.
    Sequence the arrays and returns the weights of all possible sequences!
    This version is fast but memory consuming, better for the short arrays!
    written by A. Chenani Sep. 2014
    '''
    sqStack = []
    weightsStack = []
    for item in arr:
        item.sort()
        arr_seqz = [sequencer(item)]
        #####Cunstructing sequence set
        start =  timeit.default_timer()
        
        arr_seqz_set = np.array([])
        sample_seq = arr_seqz[0].copy()
        sample_seq.sort()
        arr_seqz_set = np.append(arr_seqz_set,sample_seq)
        while next_permutation(sample_seq):
            arr_seqz_set = np.append(arr_seqz_set,sample_seq)
        arr_seqz_set = arr_seqz_set.reshape(np.math.factorial(sample_seq.size),arr_seqz[0].size)
        
        stop = timeit.default_timer()
        print 'part one --> %f' %(stop - start)
        
        
        #####creating MUA permutation list!
        start =  timeit.default_timer()
        
        while next_permutation(item):
            arr_seqz.append(sequencer(item))
        arr_seqz = np.array(arr_seqz) #This contains all sequences coming from permutations of an specific MUA!
       
        stop = timeit.default_timer()
        print 'part two --> %f' %(stop - start)
        
        #####Calculating the wights of all possible sequences
        start =  timeit.default_timer()
        weights = []
        for item in arr_seqz_set:
            selection = arr_seqz[np.where(arr_seqz[:,0]==item[0])[0]]
            for ii in range(1,arr_seqz_set.shape[1]):
                selection = selection[np.where(selection[:,ii]==item[ii])[0]]
            weights.append(selection.size/float(arr_seqz.size))
        weightsStack.append(np.array(weights))
        sqStack.append(arr_seqz_set)
        
        stop = timeit.default_timer()
        print 'part three --> %f' %(stop - start)
    return sqStack,weightsStack
def All_permutations(arr):
    '''
    produces all possible permutations of a given array using lexographical ordering.
    Sequence the arrays and returns the weights of all possible sequences!
    This version is not using that much memory but its slow, better for long arrays.
    written by A. Chenani Sep. 2014
    '''
    sqStack = []
    weightsStack = []
    for item in arr:
        item.sort()
        #####Cunstructing sequence set
        arr_seqz_set = np.array([])
        sample_seq = sequencer(item)
        sample_seq.sort()
        arr_seqz_set = np.append(arr_seqz_set,sample_seq)
        while next_permutation(sample_seq):
            arr_seqz_set = np.append(arr_seqz_set,sample_seq)
        arr_seqz_set = arr_seqz_set.reshape(np.math.factorial(sample_seq.size),sample_seq.size)
        
        
        
        
        #####creating MUA permutation list!
        start =  timeit.default_timer()
        lexRank = [0]
        while next_permutation(item):
            selection = np.arange(arr_seqz_set.size)
            for ii in range(arr_seqz_set[0].size -1 ):
                selection = np.intersect1d(selection,np.where(arr_seqz_set[:,ii]==sequencer(item)[ii])[0])
            lexRank.append(selection[0])
        lexRank = np.array(lexRank)
        stop = timeit.default_timer()
        print 'part two --> %f' %(stop - start)
        #####Calculating the wights of all possible sequences
        weights = []
        for ii,item in enumerate(arr_seqz_set):
            weights.append(np.where(lexRank == ii)[0].size /float(lexRank.size))
        weightsStack.append(np.array(weights))
        sqStack.append(arr_seqz_set)
    return sqStack,weightsStack

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
        supplied root directory.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield [path,filename]


# In[ ]:
start = timeit.default_timer()
n=6
A = range(n)
corr1T = np.array([])
corr2T = np.array([])
perms = []
tempA = range(n)
tempB = np.random.permutation(tempA)
tempA = np.random.permutation(tempA)
for item in all_permutations([range(n)])[0][0]:
    perms.append(item)
    c1 = scipy.stats.pearsonr(item,tempA)[0]
    c2 = scipy.stats.pearsonr(item,tempB)[0]
    LRcorr = np.array([c1,c2])
    corr2T = np.append(corr2T,LRcorr[np.argmax(np.abs(LRcorr))])
    corr1T = np.append(corr1T,c1)
perms = np.int8(np.array(perms))
#print perms
pkl.dump(corr2T,open('./corr'+str(n)+'.crr2','wb'),pkl.HIGHEST_PROTOCOL)
pkl.dump(corr2T,open('./corr'+str(n)+'.crr','wb'),pkl.HIGHEST_PROTOCOL)
pkl.dump(perms,open('./perms'+str(n)+'.pr','wb'),pkl.HIGHEST_PROTOCOL)
stop = timeit.default_timer()
print stop - start


