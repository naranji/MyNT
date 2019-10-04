import sys
sys.path.append('/home/chenani/ownCloud/Workspaces/Eclipse/dataAnalysis/Sleep-current/src/')
import cPickle as pkl
import pandas as pd

lfpFname = sys.argv[1]
lfp = pkl.load(open(lfpFname,'rb'))
slp = pd.read_pickle(sys.argv[2])
sessionName = lfp.tags['file'].split('/')[-2].split('-')[-1].split('p')[0]+'p0'+lfp.tags['file'].split('/')[-2].split('-')[-1].split('p')[1]
slpCurrent = slp[slp.session==sessionName]
swsEP = slpCurrent[slpCurrent.epoch=='SWS']
lfp.sws_episodes = swsEP[['t0','t1']].as_matrix()*1e3
lfp.hilbertTransform()
lfp.SWS_signal(filter=True,f_1=100,f_2=250)
lfp.ripple_recorder(rippleCut=False,rippleMix=False,SWRmix=False)
print 'Dumping %s' %lfpFname
pkl.dump(lfp,open(sys.argv[1],'wb'),pkl.HIGHEST_PROTOCOL)
print 'D>O>N>E!!!'