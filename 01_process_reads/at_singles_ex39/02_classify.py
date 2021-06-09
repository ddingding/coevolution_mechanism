
from mutTools import rev_complement, fasta_iter
from configx39 import expToInd, expToSeq
import mapClassPe as mcp


#fastaInDir1 = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/03fastqFilteredMi1/'
#dout1 = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/04_class_mi1/'
#fastaInDir2 ='/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/03fastqFilteredMi2/'
#dout2 = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/04_class_mi2/'

def classifySamples(fastaInDir, dout):
    print 'classifying...'+fastaInDir
    c=0
    for k,ind in expToInd.iteritems():
        ind= rev_complement(ind)
        template = expToSeq[k]
        with open(dout+ind+'_class.csv','w') as fout:
            for n,s in fasta_iter(fastaInDir + ind+'_merged.fasta'):
                writeL = mcp.mapAndClassify300Read(n,s, template)
                fout.write('\t'.join(writeL)+'\n')
                c+=1
                print c
        if c%10000 == 1:
            print c

fastaInDir1 = './mi1/'
dout1='./class1/'
fastaInDir2 = './mi2/'
dout2='./class2/'

classifySamples(fastaInDir1, dout1)
classifySamples(fastaInDir2, dout2)
