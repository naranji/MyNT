'''
Created on Jun 10, 2014

@author: chenani
'''
import matplotlib
matplotlib.use('Agg')
import os,sys,fnmatch
import Recordings
#from sws_rem_compare import *
import pickle as pkl
import collections
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
##########################################
animalPath = "/home/chenani/atlas-data/7714/real"
##########################################FUNCTIONS
def locate(pattern, root=os.curdir):
    '''
    Locate all files matching supplied filename pattern in and below
        supplied root directory.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield [path,filename]
def save(report, filename):
        '''
        This function saves the ephys object to a .ephys file!
        '''
        pkl.dump(report, open(filename + ".rpt", 'wb'), pkl.HIGHEST_PROTOCOL)
        print "data has been saved to %s.rpt" % filename
##########################################
ephysList = []
for ephys in locate("*.ephys",animalPath):
    ephysList.append(ephys)
#########################################
for item in ephysList:
    print os.path.join(item[0],item[1])
    ephysFile = os.path.join(item[0],item[1])
    dataFolder = item[0]+"/catalog"
    corrupted = 0.0
    corr_percent =  []
    abs_sig_avg = []
    fileNames = []
    qualityRank = []
    strengthRank = []
    if os.path.isfile(os.path.join(item[0],'catalog/report.rpt')):
        print "The report file alarm.."
    rec = pkl.load(open(ephysFile, "rb"))
    ######################################## Operations on different channels, picking the "best" csc!
    for item in rec.LFPs:
        corrupted += len(np.where(item.signal > 32766 )[0])
        corrupted += len(np.where(item.signal < -32766 )[0])
        corr_percent.append(round(corrupted / item.signal.size * 100,3))
        abs_sig_avg.append(np.abs(item.signal).mean())
        fileNames.append(item.tags['file'])
    sig_length = round(rec.LFPs[0].signal.size * rec.LFPs[0].dt / 60000,1)
    corr_percent = np.array(corr_percent)
    abs_sig_avg = np.array(abs_sig_avg)
    for i in range(corr_percent.size):
        qualityRank.append(np.where(corr_percent < corr_percent[i])[0].size + 1)
        strengthRank.append(np.where(abs_sig_avg > abs_sig_avg[i])[0].size + 1)
    report = collections.OrderedDict()
    report['name'] = fileNames
    report['duration (min)'] =  sig_length
    report['quality_rank'] = np.array(qualityRank)
    report['strength_rank'] =  np.array(strengthRank)
    csc_rank = report['quality_rank']+report['strength_rank']
    recommended_csc = fileNames[csc_rank.argmin()]
    csc = rec.LFPs[csc_rank.argmin()]
    ##################################### Generating numbers and plots for the report.
    if  not os.path.exists(dataFolder):
        os.mkdir(dataFolder)
    ######
    if not hasattr(csc,'signal_filtered'): #and (csc.cutfreqz[0] < 140 or csc.cutfreqz[0] == 250 ):
        csc.filter(140,250)
    if not hasattr(csc,'sws_signal'):
        csc.SWS_signal(True)
    if not hasattr(csc,'rem_signal'):
        csc.REM_signal()
    if  not hasattr(csc,'hilbertAbsolute'):
        csc.hilbertTransform()
    if  not hasattr(csc,'ripples'):
        csc.ripple_recorder(removeREMripples=True)
    if not hasattr(csc,'IRI'):
        csc.iri()
    
    ############Removing ripples during REM
    ripp_copy = csc.ripples
    SUM = 0
    for item in csc.rem_episodes:
        j = 0
        for jtem in csc.ripples:
            peak_in = jtem[-1] <item[1] and jtem[-1]>item[0]
            starts_in = jtem[0] <item[1] and jtem[0]>item[0]
            ends_in = jtem[1] <item[1] and jtem[1]>item[0]
            if peak_in or starts_in or ends_in :
                ripp_copy = np.delete(ripp_copy,j-SUM,0)
                SUM += 1
                j+= 1
    if csc.rem_signal.shape[0] > 1:
        REM_percentile = round(np.hstack(csc.rem_signal).size * csc.dt / csc.duration() *100,1)
        REM_ripp = round(SUM / float(csc.ripples.shape[0]) *100,1)
    else:
        REM_percentile = 0
        REM_ripp = 0
    ####################################PLOTS!
    with PdfPages(dataFolder + 'plots.pdf') as pdf:
        f1 = pl.figure()
        ax1 = f1.add_subplot(111)
        ripp_duration = dataFolder+'/ripp_duration.png'
        csc.ripple_statistics(fig=f1,ax1=ax1)
        f1.savefig(ripp_duration)
        pdf.savefig()
        pl.close()
        f2 = pl.figure()
        ax2 = f2.add_subplot(111)
        ripp_iri = dataFolder+'/iri.png'
        csc.ripple_statistics(fig=f2,ax3=ax2)
        f2.savefig(ripp_iri)
        pdf.savefig()
        pl.close()
        f3 = pl.figure()
        ax3 = f3.add_subplot(111)
        ripp_fft = dataFolder+'/ripp_fft.png'
        csc.ripple_statistics(fig=f3,ax4=ax3)
        f3.savefig(ripp_fft)
        pdf.savefig()
        pl.close(f3)
        fig = pl.figure()
        ax = fig.add_subplot(111)
        rec.LFPs.plot(fig,ax)
        sig_path = dataFolder+'/signal.png'
        fig.savefig(sig_path)
        pl.close(fig)
        
    
    print "========================================"
    report['recommended_csc'] =  recommended_csc
    report['corrupted_signal_percentile'] = corr_percent
    report['rem_percentile'] = REM_percentile
    report['rem_ripples_percentile'] =  REM_ripp
    report['number of ripples'] = ripp_copy.shape[0]
    report['IRI_mean (ms)'] = csc.IRI.mean()
    report['IRI_std (ms)'] = csc.IRI.std() 
    report['ripples'] = ripp_copy
    
    keys = report.keys()
    with open(dataFolder+"/report.txt", "w") as f:
        for item in keys:
            print item
            f.write('================='+item + '\n')
            line = str(report.get(item)) + '\n'
            f.write(line)
        f.close()
    save(report,dataFolder+'/report')
    rec.save(ephysFile.split('.')[0])
    del rec