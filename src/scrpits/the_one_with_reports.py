'''
Created on Jun 11, 2014

@author: chenani
'''

import os,sys,fnmatch
import pickle as pkl
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages

##########################################
animalPath = "/home/chenani/atlas-data/ali/8981"
##########################################FUNCTIONS
histograms = False
barPlots = False
writeText = False
for argv in sys.argv[1:]:
    if argv == '-hist':
        histograms = True
    if argv == '-bars':
        barPlots = True
        
def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
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
########################################## Creating report list
rptList = []
for rpt in locate("*.rpt",animalPath):
    rptList.append(rpt)
######################################### Writing reports to .txt files.
if writeText:
    for item in rptList:
        if str(item).find('leep') > 1:
            print os.path.join(item[0],item[1])
            rptFile = os.path.join(item[0],item[1])
            report = pkl.load(open(rptFile, "rb"))
            #swr_rate =  report['number of ripples']/report['duration']/60
#             report['swr_mean_rate'] = swr_rate
            keys = report.keys()
            with open(rptFile.split('.')[0]+".txt", "w") as f:
                for item in keys:
                    print item
                    f.write('================='+item + '\n')
                    line = str(report.get(item)) + '\n'
                    f.write(line)
                f.close()
            save(report,rptFile.split('.')[0])
###############################################
vrstack = np.array([])
ltstack = np.array([])
compstack = np.array([])
acceptableDates = np.array(['2014-06-16','2014-07-28','2014-07-29','2014-07-30',
                   '2014-08-04','2014-08-05','2014-08-06','2014-08-07'])
if histograms:
    for item in rptList:
        if str(item).find('9588/vr') > 1 and str(item).find('leep2') > 1:
            print os.path.join(item[0],item[1])
            rptFile = os.path.join(item[0],item[1])
            report = pkl.load(open(rptFile, "rb"))
            vrstack = np.append(vrstack,report['ripples'])
            vrstack = vrstack.reshape(vrstack.shape[0]/4,4)
        elif str(item).find('9588/real') > 1 and str(item).find('leep2') > 1:
            print os.path.join(item[0],item[1])
            rptFile = os.path.join(item[0],item[1])
            report = pkl.load(open(rptFile, "rb"))
            ltstack = np.append(ltstack,report['ripples'])
            ltstack = ltstack.reshape(ltstack.shape[0]/4,4)
            compstack = np.append(compstack,report['ripples'])
        
        elif str(item).find('Sleep_after_realTrack') > 1 and str(item).find('leep') > 1:
            print os.path.join(item[0],item[1])
            rptFile = os.path.join(item[0],item[1])
            report = pkl.load(open(rptFile, "rb"))
            compstack = np.append(compstack,report['ripples'])
            compstack = compstack.reshape(compstack.shape[0]/4,4)
        vr_duration = vrstack[:,1] - vrstack[:,0]
        vr_iri = np.diff(vrstack[:,3]) 
        vr_iri = vr_iri[np.where(vr_iri < 5000)[0]]
        lt_iri = np.diff(ltstack[:,3])
        lt_iri = lt_iri[np.where(np.abs(lt_iri) < 5000)[0]];
        lt_duration = ltstack[:,1] - ltstack[:,0]
        comp_du = compstack[:,1] - compstack[:,0]
        comp_ir = np.diff(compstack[:,3])
        comp_ir = comp_ir[np.where(np.abs(comp_ir) < 5000)[0]];
if barPlots:
    rptList2 = []
    for pattern in acceptableDates:
        for item in rptList:
            if str(item).find(pattern) > 1:
                rptList2.append(item)
    
    for item in rptList:
        if str(item).find('real') > 1 and str(item).find('leep2') > 1:
            print os.path.join(item[0],item[1])
            rptFile = os.path.join(item[0],item[1])
            report = pkl.load(open(rptFile, "rb"))
            vrstack = np.append(vrstack,report['ripples'])
            vrstack = vrstack.reshape(vrstack.shape[0]/4,4)
        elif str(item).find('vr') > 1 and str(item).find('leep2') > 1:
            print os.path.join(item[0],item[1])
            rptFile = os.path.join(item[0],item[1])
            report = pkl.load(open(rptFile, "rb"))
            ltstack = np.append(ltstack,report['ripples'])
            ltstack = ltstack.reshape(ltstack.shape[0]/4,4)
            compstack = np.append(compstack,report['ripples'])
    lt_duration = ltstack[:,1] - ltstack[:,0]
    vr_duration = vrstack[:,1] - vrstack[:,0]
    lt_iri = np.diff(ltstack[:,3])
    lt_iri = lt_iri[np.where(lt_iri > 0)[0]]
    lt_iri = lt_iri[np.where(lt_iri < 5000)[0]]
    vr_iri = np.diff(vrstack[:,3])
    vr_iri = vr_iri[np.where(vr_iri > 0)[0]]
    vr_iri = vr_iri[np.where(vr_iri < 5000)[0]]
    #############################################PLOTS
    fig = pl.figure()
    ax = fig.add_subplot(111)

    ## the data
    N = 2
    lt_mean = [lt_iri.mean(),lt_duration.mean()]
    lt_std =  [lt_iri.std(),lt_duration.std()]
    vr_mean = [vr_iri.mean(),vr_duration.mean()]
    vr_std =  [vr_iri.std(),vr_duration.std()]

    ## necessary variables
    ind = np.arange(N)                # the x locations for the groups
    width = 0.35                      # the width of the bars

    ## the bars
    rects1 = ax.bar(ind, lt_mean, width,
                    color='green',
                    yerr=lt_std,
                    error_kw=dict(elinewidth=2,ecolor='gray'))

    rects2 = ax.bar(ind+width, vr_mean, width,
                    color='blue',
                    yerr=vr_std,
                    error_kw=dict(elinewidth=2,ecolor='gray'))
    xTickMarks = ['IRI', 'Duration']
    xtickNames = ax.set_xticklabels(xTickMarks)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind+width)
    pl.setp(xtickNames, rotation=45, fontsize=10)

    ## add a legend
    ax.legend( (rects1[0], rects2[0]), ('REAL', 'VR') )

    pl.show()
    

jj = 0
print '------------------------------------'
du_sdev = np.array([])
ir_sdev = np.array([])
du_mean = np.array([])
ir_mean = np.array([])

# while jj < comp_du.size / 100 :
#     du_mean = np.append(du_mean,comp_du[range(100* jj,100* jj+100)].mean())
#     du_sdev = np.append(du_sdev,comp_du[range(100* jj,100* jj+100)].std())
#     jj+=1
# jj = 0
# while jj < comp_ir.size / 100 :
#     ir_mean = np.append(ir_mean,comp_ir[range(100* jj,100* jj+100)].mean())
#     ir_sdev = np.append(ir_sdev,comp_ir[range(100* jj,100* jj+100)].std())
#     jj+=1
#  
# print ir_mean
# print ir_sdev
# with PdfPages('//home/chenani/Pictures/Plots/histograms.pdf') as pdf:
#     duration_fig =pl.figure()
#     duration_fig.suptitle('SWR duration')
#     vr_ax_du = duration_fig.add_subplot(111)
#     lt_ax_du = duration_fig.add_subplot(111)
#     iri_fig =pl.figure()
#     iri_fig.suptitle('Inter SWR interval')
#     vr_ax_iri = iri_fig.add_subplot(111)
#     lt_ax_iri = iri_fig.add_subplot(111)
#     delx_du = 2 *du_mean.mean()
#     vr_ax_du.hist(vr_duration, bins=100,alpha=0.95,histtype='stepfilled',label='VR',log=True)
#     lt_ax_du.hist(lt_duration,bins=100,alpha=0.95,histtype='stepfilled',label='real',log=True)
#     lt_ax_du.axvspan(3*(du_mean.mean()+du_mean.std()) - delx_du,3*(du_mean.mean()-du_mean.std()) - delx_du,0.8,0.9, color='r', linewidth=1,alpha = 0.4)
#     vr_ax_du.axvline(3*vr_duration.mean() - delx_du,0.82,0.88, color='b', linewidth=3)
#     lt_ax_du.axvline(3*lt_duration.mean() - delx_du,0.82,0.88, color='g', linewidth=3)
#     lt_ax_du.axvline(du_mean.mean(),0.8,0.9, color='k', linewidth=1)
#     lt_ax_du.set_xlabel('Time(ms)')
#     #lt_ax_du.set_xlim(30,200)
#     lt_ax_du.legend(loc='upper right')
#     delx_ir = 2 *ir_mean.mean()
#     vr_ax_iri.hist(vr_iri, bins=100,alpha=0.95,histtype='stepfilled',label='VR',log=True)
#     lt_ax_iri.hist(lt_iri,bins=100,alpha=0.95,histtype='stepfilled',label='real',log=True)
#     lt_ax_iri.axvline(3*ir_mean.mean(),0.82,0.88, color='g', linewidth=2)
#     lt_ax_iri.axvspan(3*(ir_mean.mean()+ir_mean.std()) - delx_ir,3*(ir_mean.mean()-ir_mean.std()) - delx_ir,0.8,0.9, color='r', linewidth=1,alpha = 0.4)
#     vr_ax_iri.axvline(3*vr_iri.mean() - delx_ir,0.82,0.88, color='b', linewidth=3)
#     lt_ax_iri.axvline(3*lt_iri.mean() - delx_ir,0.82,0.88, color='g', linewidth=3)
#     lt_ax_iri.axvline(ir_mean.mean(),0.8,0.9, color='k', linewidth=1)
#     lt_ax_iri.set_xlabel('Time(ms)')
#     #pl.xlim(0,2000)
#     lt_ax_iri.legend(loc='upper right')
#     pdf.savefig()
#     pl.show()
#     print vr_duration.size
