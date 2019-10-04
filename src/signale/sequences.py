"""
Signale.sequences
A module for sequence related calculations
"""
__author__ = ("A. Chenani")
__version__ = "1.0, May 2014"
# python modules

# other modules
import numpy as np
from itertools import combinations
# custom made modules

# package modules

###################################################### FUNCTIONS
def p1(j,m):
    """
    This function computes the probability of observing a sequence by chance, having a finite set of neurons with
    different firing rates. For more details consult Smith et. al. 2010...
  
    Parameters
    j: length of the sequnce. (int)
    m: firing rate vector, consisted of firing rate of k'th unit as it's kth element.
  
    Output:
    It gives you a floating point number indicating the probability!!!
    """
    N = len(m)  # Number of Units
    l = np.array(list(combinations(range(N),j))) # array of all possible combinations of length j
    nt = np.sum(m) 
    sum = 0
    for item in l:
        sum += np.prod(m[item])
    return float(sum)/nt**j    
###################################################### CLASSES
