"""
signale.tools
=============

A module for general purpose tools.
"""

__author__ = ("KT")
__version__ = "1.2, August 2013"

# python modules

# other modules
import numpy
import os,fnmatch
# custom made modules

# package modules




###################################################### FUNCTIONS


def findMaximaMinima(a, findMaxima=1, findMinima=1, or_equal=1):
    """ Find the maxima and minima in an 1D array.

    Paramters:
    a ... 1D array
    findMaxima, findMinima ... boolean, if ought to look for maxima and/or minima [defaults 1]

    Returns:
    dictionary{'maxima_number', 'minima_number', 'maxima_indices', 'minima_indices'}
    """
    gradients = numpy.diff(a)

    maxima_num = 0
    minima_num = 0
    max_indices = []
    min_indices = []

    if or_equal:  # includes only the next entries and they can also be equal to zero
        for i, g in enumerate(gradients[:-1]):
            if findMaxima and cmp(g, 0) >= 0 and cmp(gradients[i + 1], 0) <= 0 and g != gradients[i + 1]:
                maxima_num += 1
                max_indices.append(i + 1)

            if findMinima and cmp(g, 0) <= 0 and cmp(gradients[i + 1], 0) >= 0 and g != gradients[i + 1]:
                minima_num += 1
                min_indices.append(i + 1)
    else:  # includes also the previous and next entries and also expects them not to be equal to zero
        for i, g in enumerate(gradients[1:-1]):
            j = i + 1  # get real index
            if findMaxima and cmp(g, 0) >= 0 and cmp(gradients[j - 1], 0) > 0 and cmp(gradients[j + 1], 0) < 0:
                maxima_num += 1
                max_indices.append(j + 1)

            if findMinima and cmp(g, 0) <= 0 and cmp(gradients[j - 1], 0) < 0 and cmp(gradients[j + 1], 0) > 0:
                minima_num += 1
                min_indices.append(j + 1)

    return {'maxima_number': maxima_num , 'minima_number': minima_num, \
            'maxima_indices': numpy.array(max_indices), 'minima_indices': numpy.array(min_indices)}


def findNearest(array, value):
    """ Find the entry nearest to a value in an 1D array.

    Paramters:
    array ... 1D array
    value ... value to find nearest entry for

    Returns:
    a tuple (index, array[index])
    """
    index = (numpy.abs(array - value)).argmin()
    return (index, array[index])


def findNearest_vec(a, v):
    """ Find the entry nearest to a vector value in an 1D array.

    Paramters:
    a ... 1D array
    v ... value to find nearest entry for

    Returns:
    a tuple (index, a[index])
    """
    diff = a - v
    norms = numpy.array([])
    for d in diff:
        norms = numpy.append(norms, numpy.linalg.norm(d))
    index = norms.argmin()
    return (index, a[index])



def sameLengthArr(arr):
    ''' Reduce a list with several 1D arrays to the length of the smallest.

        Expects a list with several 1D arrays. The
        arrays may be of different size.

        parameters:
        arr ... list of 1D arrays

        Returns:
        numpy 2D Array
    '''

    minLength = 10 ** 100
    for a in arr:
        minLength = min(minLength, a.size)

    return numpy.array([a[:minLength] for a in arr])


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len < 3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s = numpy.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]  # enlarge signal at the borders
    #    print x.shape, s.shape

    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='same')
    y = y[window_len - 1:-(window_len - 1)]  # return to original signal size
    #   print y.shape

    return y


def transposeArr(arr, numeric=True):
    ''' Transpose a list of 1D arrays.

        Expects a list with several 1D arrays. The
        arrays may be of different size.

        parameters:
        arr ... list of 1D arrays
        numeric ... True if the indidivual arrays (lists) contain numbers,
                    then a list of numpy arrays is returned

        Returns:
        transposed list
    '''

    maxLength = 0
    for a in arr:
        maxLength = max(maxLength, a.size)

    transposedArr = []
    for i in range(maxLength):
        l = []
        for a in arr:
            if i < a.size:
                l.append(a[i])
        if numeric:
            transposedArr.append(numpy.array(l))
        else:
            transposedArr.append(l)

    return transposedArr

def locate(pattern, root=os.curdir):
    '''
    Locate all files matching supplied filename pattern in and below
        supplied root directory.
    Parameters
    ----------
        pattern : string
            A string representing pattern you want to look for in the path.
        root : string
            contains the path address you want to look trough!
    Returns
    ----------
        Array of strings consist of two column, first contains the path and second the file names that matched
        the desired pattern.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield [path,filename]


###################################################### TOOLS



###################################################### CLASSES
