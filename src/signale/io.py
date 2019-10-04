"""
signale.io
==========

A module for input/output functions for signals. 
"""
__author__ = ("KT")
__version__ = "4.0, May 2013"


# python modules
import struct, os, sys

# other modules
import numpy

# custom made modules

# package modules
import cscs
import place_cells
import signals
import timeit as tit

###################################################### FUNCTIONS


def load_tFile(fileName, showHeader=False):
    """
    For reading AD Redishs's MClust 3.5 t-files with python/numpy.
    These files contain events/spikes data.
    
    Returns a Neurotools.signals SpikeTrain object. 
    """
    spikes = []
    header = ''
    with open(fileName, 'rb') as f:
        headerFinished = False
        while not headerFinished:
            line = f.readline()
            if line.startswith('%%ENDHEADER'):
                headerFinished = True
            header += line
        while True:
            chunk = f.read(4)
            if chunk:
                spike = struct.unpack('>L', chunk)[0]  # data in MClust t files is stored as UInt32 big endian
                spikes.append(spike)
            else:
                break
    spikes = place_cells.placeCell_spikezug(numpy.float_(spikes) / 10.)  # convert to ms, 
                                                               # time in MClust t-files is in tens of ms
    
    if showHeader:
        print header

    return spikes

def load_ncsFile(fileName, showHeader=False):
    """
    For loading Neuralynx continuously sampled channel (CSC) recorded data
    with python/numpy.
    Returns a NeuralynxCSC object.
    
    Notes:
    ----------
        * It's faster than load_ncs_File_old by a factor of 50!
        *There is a periodic difference of 0.06 ms with the 
        timeAxis calculated from _calcTimeAxis method with the period 
        of NumvalidSamples.
        :rtype: NeuralynxCSC

        
    """
    with open(fileName, 'rb') as f:
    
        header = f.read(16 * 2 ** 10)  # the first 16 kB are an ASCII text header
    
    if showHeader:
        print header
    
    ncs_dtype = numpy.dtype([    
        ('TimeStamps'  ,     '<u8'),
        ('ChannelNumber'  , '<u4'),
        ('SampleFreq',      '<u4'),
        ('NumValidSamples', '<u4'),
        ('Samples'  ,       '<i2',(512,)),
    ])
    memmap = numpy.memmap(fileName,dtype=ncs_dtype,offset=(16 * 2**10)) 
    
    timeStamps = numpy.array(memmap[:]['TimeStamps'])  # change to ms, Neuralynx time stamps are in microseconds
    dataSamples = numpy.array(memmap[:]['Samples']).flatten()
    sampleFreq = memmap[0]['SampleFreq']
    timeAxis = numpy.array([])
    dt = 1000.0/sampleFreq

     
    csc = cscs.NeuralynxCSC(timeStamps, dataSamples, sampleFreq)
    csc.tags['file'] = fileName
    csc.tags['path'] = os.getcwd()
    csc.tags['channel'] = memmap[0]['ChannelNumber']
    csc.tags['validSamples'] = numpy.array(memmap[:]['NumValidSamples'])
    csc.header = header
    # a = header.find('ADBitVolts')
    # csc.ADBitVolts = float(header[a + 11:a +37])      #AD Unit in Volts! 
    csc.dataSamplesNums = numpy.array([memmap[item]['Samples'].size for item in range(len(memmap))]) #####
    
    ti = tit.time.time()
    if not numpy.diff(numpy.array(memmap[:]['NumValidSamples'])).sum():
        Nvalid = memmap[0]['NumValidSamples']
        tChunk = numpy.arange(0,Nvalid*1000*dt,1000*dt)
        tBlock = numpy.tile(tChunk,(timeStamps.size,1)) 
        tBlock = numpy.array([tBlock[ii,:]+ timeStamps[ii] for ii in range(timeStamps.size)])
        timeAxis = tBlock.flatten() / 1000.
    else:
        for ii in range(timeStamps.size):
            timeAxis= numpy.append(timeAxis,numpy.arange(0,memmap[ii]['NumValidSamples']*1000*dt,1000*dt)+timeStamps[ii])
    print 'TimeAxis calculation took %s seconds' %(tit.time.time() - ti)
    
    csc.timeAxis = timeAxis
    if numpy.any(numpy.diff(numpy.diff(timeStamps))):
        print('WARNING: Not all records are contiguous. Packets lost?')
    if (csc.dataSamplesNums - csc.tags['validSamples']).sum():
        print 'WARNING: Number of samples and number of valid samples not the same.'
    return csc   
    
    

def load_ncsFile_old(fileName, showHeader=False):
    """
    For loading Neuralynx continuously sampled channel (CSC) recorded data
    with python/numpy.
    Returns a NeuralynxCSC object. 
    """
    timeStamps = []
    dataSamples = []
    dataSamplesNums = []
    validSamples = []
    with open(fileName, 'rb') as f:
    
        header = f.read(16 * 2 ** 10)  # the first 16 kB are an ASCII text header
        
        if showHeader:
            print header

        count = -1
        while True:
            qwTimeStamp = f.read(8)  # UInt64, Cheetah timestamp for this record. This value is in microseconds.
            dwChannelNumber = f.read(4)  # UInt32, channel number for this record.
            dwSampleFreq = f.read(4)  # UInt32, sampling frequency reported by the acquisition hardware when recording this record.
            dwNumValidSamples = f.read(4)  # UInt32, Number of values in snSamples containing valid data.
            snSamples = f.read(512 * 2)  # Int16[], Data points for this record.

            if qwTimeStamp:
                qwTimeStamp = struct.unpack('L', qwTimeStamp)[0]
                dwChannelNumber = struct.unpack('I', dwChannelNumber)[0]
                dwSampleFreq = struct.unpack('I', dwSampleFreq)[0]
                dwNumValidSamples = struct.unpack('I', dwNumValidSamples)[0]
                
                channel = dwChannelNumber
                sampleFreq = dwSampleFreq
                validSamples.append(dwNumValidSamples)
                
                snSamples = [snSamples[i:i + 2] for i in range(0, snSamples.__len__(), 2)]
                samples = []
                for sample in snSamples:
                    sample = struct.unpack('h', sample)[0]
                    samples.append(sample)
                snSamples = numpy.array(samples)
                
                timeStamps.append(qwTimeStamp)
                dataSamples.append(snSamples)
                dataSamplesNums.append(snSamples.size)
                
                if dataSamplesNums[-1] != validSamples[-1]:
                    print 'WARNING: Number of samples and number of valid samples not the same.'
                
                
                # print some data?
                if count < 0 and showHeader:
                    count += 1
                    print '----> count', count      
                    print qwTimeStamp         
                    print channel          
                    print sampleFreq             
                    print validSamples[-1]                   
                    # print snSamples
                    print snSamples.shape              
                    print ''
            else:
                break
    
    timeStamps = numpy.array(timeStamps) / 1000.  # change to ms, Neuralynx time stamps are in microseconds
    dataSamples = numpy.array(dataSamples)
    dataSamples = dataSamples.flatten()
    
    # integrity check, after extracellpy by Santiago Jaramillo
    if numpy.any(numpy.diff(numpy.diff(timeStamps))):
        print('WARNING: Not all records are contiguous. Packets lost?')
    
    csc = cscs.NeuralynxCSC(timeStamps, dataSamples, sampleFreq)
    csc.tags['file'] = fileName
    csc.tags['path'] = os.getcwd()
    csc.tags['channel'] = channel
    csc.tags['validSamples'] = numpy.array(validSamples)
    csc.header = header
    # a = header.find('ADBitVolts')
    # csc.ADBitVolts = float(header[a + 11:a +37])      #AD Unit in Volts! 
    csc.dataSamplesNums = dataSamplesNums
    ti = tit.time.time()
    csc.recalc_timeAxis()  # recalculate time values and timeAxis using the validSamples
    print 'TimeAxis calculation took %s seconds' %(tit.time.time() - ti)

    return csc


def load_nevFile(fileName, showHeader=False):
    """
    For loading Neuralynx Event files with python/numpy.
    Returns a NeuralynxEvent object. 
    """
    timeStamps = []
    eventStrings = []
    
    with open(fileName, 'rb') as f:
    
        header = f.read(16 * 2 ** 10)  # the first 16 kB are an ASCII text header
        
        if showHeader:
            print header
        
        count = 0
        while True:
            nstx = f.read(2)  # Int16, Reserved
            npkt_id = f.read(2)  # Int16, ID for the originating system of this packet.
            npkt_data_size = f.read(2)  # Int16, This value should always be two (2).
            qwTimeStamp = f.read(8)  # UInt64, Cheetah timestamp for this record. This value is in microseconds.
            nevent_id = f.read(2)  # Int16, ID value for this event.
            nttl = f.read(2)  # Int16, Decimal TTL value read from the TTL input port.
            ncrc = f.read(2)  # Int16, Record CRC check from Cheetah. Not used in consumer applications.
            ndummy1 = f.read(2)  # Int16, Reserved
            ndummy2 = f.read(2)  # Int16, Reserved
            dnExtra = f.read(8 * 4)  # Int32[], Extra bit values for this event. 
                                                # This array has a fixed length of eight (8).
            eventString = f.read(128)  # Event string associated with this event record. This string has a maximum length of 128 characters.

            if nstx:
                nstx = struct.unpack('h', nstx)[0]
                npkt_id = struct.unpack('h', npkt_id)[0]
                npkt_data_size = struct.unpack('h', npkt_data_size)[0]
                qwTimeStamp = struct.unpack('L', qwTimeStamp)[0]
                nevent_id = struct.unpack('h', nevent_id)[0]
                nttl = struct.unpack('h', nttl)[0]
                ncrc = struct.unpack('h', ncrc)[0]
                ndummy1 = struct.unpack('h', ndummy1)[0]
                ndummy2 = struct.unpack('h', ndummy2)[0]
                
                timeStamps.append(qwTimeStamp)
                eventStrings.append(eventString.split('\x00')[0])
                
                # print some data?
                if count < 0:
                    count += 1
                    print '----> count', count 
                    print nstx                    
                    print npkt_id                
                    print npkt_data_size        
                    print qwTimeStamp         
                    print nevent_id          
                    print nttl             
                    print ncrc                   
                    print ndummy1              
                    print ndummy2                
                    print dnExtra              
                    print eventString
                    print ''
            else:
                break
        
        timeStamps = numpy.array(timeStamps) / 1000.  # change to ms, Neuralynx time stamps are in microseconds
        
        nev = signals.NeuralynxEvents(timeStamps, eventStrings)
        nev.tags['file'] = fileName
        nev.tags['path'] = os.getcwd()
        
        return nev



def load_rawFile(fileName, exclude=[], showHeader=False):
    """
    For loading RAW data file exported from Multi Channel Systems MCRack with python/numpy.
    """
    data = []
    dataSamplesNums = []
    
    header = ''
    with open(fileName, 'rb') as f:
        headerFinished = False
        while not headerFinished:
            line = f.readline()
            if line.startswith('EOH'):
                headerFinished = True
            header += line
            
            # get metadata from header
            if line.startswith('Sample rate'):
                sampleFreq = float(line.split('=')[1])
            elif line.startswith('Streams'):
                channels = line.split('=')[1]
                channels = channels.split(';')
                channels[-1] = channels[-1][:-2]  # the last two entries are special characters
                numChannels = channels.__len__()
        
        while True:
            chunk = f.read(4)
            if chunk:
                dataSample = struct.unpack('I', chunk)[0]  # data as UInt16
                data.append(dataSample)
            else:
                break
    
    # split data from file into the different channels
    # the raw data is stored such that at each time step the entires for all channels come one after the other
    # i.e. t1: d_1, d_2, ..., d_n 
    #      t2: d_1, d_2,..., d_n
    #      t3:  d_1, d_2,..., d_n
    ID = -1
    dataList = cscs.CSCsignalList()
    for i, ch in enumerate(channels):
        if not i + 1 in exclude:
            d = cscs.CSCsignal(data[i::channels.__len__()], sampleFreq)
            d.tags['channel'] = ch
            d.tags['file'] = fileName
            d.tags['path'] = os.getcwd()
        
            ID += 1
            dataList.append(ID, d)
 
    
    # add metadata to tags
    dataList.tags['file'] = fileName
    dataList.tags['path'] = os.getcwd() 
    dataList.tags['numChannels'] = numChannels
    dataList.tags['channels'] = channels
    dataList.header = header

    if showHeader:
        print header
    
    return dataList


def _read_metadata(fileName, convert2variables=False, showHeader=False):
    """
    Read the informations that may be contained in the header of
    the file.
    """
    metadata = {}
    cmd = ''
    f = open(fileName, 'r')
    for line in f.readlines():
        if line[0] == '#':
            cmd += line[1:].strip() + ';'
            if showHeader:
                print line.strip()
        elif line[0] == '%' and showHeader:  # for special comments that may be included in the header
            print line.strip()
        else:
            break
    f.close()
    exec cmd in None, metadata
    
    if showHeader:
        print ''    
    
    if convert2variables:
        for key in metadata.keys():
            exec key + '=' + str(metadata[key])

    return metadata


###################################################### CLASSES

