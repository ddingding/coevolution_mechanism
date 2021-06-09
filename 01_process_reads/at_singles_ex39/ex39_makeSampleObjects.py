from configx39 import expToInd, exp_to_od, expToTemplate

from os import listdir
from os.path import isfile, join
import pickle
import libClass as lc
from mutTools import rev_complement
import plottingTools as pt


####################################make pickles of the samples
##############################################################################################
def makeSamplePickles(classPath1, classPath2, pickleDir, pickleN=''):
    for k,ind in expToInd.iteritems():
        ind=rev_complement(ind)
        if not isfile(pickleDir+ind+pickleN+'.p'):
            #self, filename,sampleName=None, timepoint=None, od=None,template=None, use_nns = False, at = None, t = None)

            sampletype, timepoint = k.split('_')
            od = exp_to_od[k]
            template = expToTemplate[k]
            #BUG used classPath1 twice!! now changed to classPath2 in the second instance
            currObj = lc.SampleObject([classPath1+ind+'_class.csv',classPath2+ind+'_class.csv'], sampleName=ind, timepoint=int(timepoint), template =template, od=od, use_nns=False)
            pickle.dump(currObj, open(pickleDir+ind+pickleN+'.p', 'wb'))
        else:
            continue

#combine sample objects from both miseqs
pickleDir = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/pickles/'
makeSamplePickles(classPath1, classPath2, pickleDir, pickleN = '_miMerged_nov2019')



#######################################################################################################


pickleDir = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/pickles/'
# load all the s_objects
indToSampleObj= {}

print 'making tc'

def makeTcObj(expList, piName, template):

    for ind in [expToInd[exp] for exp in expList]: #[expToInd['1.1_0'],expToInd['1.1_330'],expToInd['1.1_600']]:
        #load samples
        print 'loading sample', ind
        ind_rev = rev_complement(ind)
        currObj = pickle.load(open(pickleDir+ind_rev+'_miMerged_nov2019.p', 'rb'))
        indToSampleObj[ind] = currObj

    #make sampleList
    sObjList = [indToSampleObj[expToInd[exp]] for exp in expList]
    tc = lc.Timecourse(template, sObjList)
    tc.calculateFitness(first_timepoint=0)
    tc.samples[-1].calc_average_aa_fit()
    pickle.dump(tc, open(pickleDir+piName, 'wb'))


makeTcObj(['4.1_0','4.1_600'], '41_tc_nov2019.p', 'pard')
makeTcObj(['4.2_0','4.2_600'], '42_tc_nov2019.p', 'pard')
