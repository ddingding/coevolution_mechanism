# takes a file of classified reads, and create sample and timecourse objects based on metadata in config file

from os.path import isfile
import pickle

from configx39 import expToInd, exp_to_od, expToTemplate
import libClass as lc
from mutTools import rev_complement


####################################make pickles of the samples
##############################################################################################

def makeSamplePickles(classPaths, pickleDir, pickleN=''):
    # go through antitoxin conditions and respective illumina multiplex indices
    # creates sample object pickles named after the illumina multiplex index used
    for k, ind in expToInd.iteritems():
        ind = rev_complement(ind)
        pickle_fname =pickleDir + ind + pickleN + '.p'
        if not isfile(pickle_fname):
            # fetch the metadata
            sampletype, timepoint = k.split('_')
            od = exp_to_od[k]
            template = expToTemplate[k]

            # create sample object
            currObj = lc.SampleObject([fin + ind + '_class.csv' for fin in classPaths], sampleName=ind,
                                      timepoint=int(timepoint), template=template, od=od, use_nns=False)
            pickle.dump(currObj, open(pickle_fname, 'wb'))
        else:
            continue

#make sample object from files for miseq run 1
classPath1 = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/04_class_mi1/'
pickleDir1 = './pickles/miseq1/'
makeSamplePickles(classPath1, pickleDir1)

#make sample object from files for miseq run 2
classPath2 = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/04_class_mi2/'
pickleDir2 = './pickles/miseq2/'
makeSamplePickles(classPath2, pickleDir2)

# combine sample objects from both miseqs
pickleDir = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/pickles/'
makeSamplePickles(classPath1, classPath2, pickleDir, pickleN='_miMerged_nov2019')

#######################################################################################################


pickleDir = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/pickles/'

print 'making tc'
def makeTcObj(expList, piName, template):
    #make timecourse object pickles

    # load all the s_objects
    indToSampleObj = {}
    for ind in [expToInd[exp] for exp in expList]:
        # load samples
        print 'loading sample', ind
        ind_rev = rev_complement(ind)
        currObj = pickle.load(open(pickleDir + ind_rev + '_miMerged_nov2019.p', 'rb'))
        indToSampleObj[ind] = currObj

    # fetch sample objects belonging to one timecourse experiment, and create a timecourse object
    sObjList = [indToSampleObj[expToInd[exp]] for exp in expList]
    tc = lc.Timecourse(template, sObjList)
    pickle.dump(tc, open(pickleDir + piName, 'wb'))


makeTcObj(['4.1_0', '4.1_600'], '41_tc_nov2019.p', 'pard')
makeTcObj(['4.2_0', '4.2_600'], '42_tc_nov2019.p', 'pard')
