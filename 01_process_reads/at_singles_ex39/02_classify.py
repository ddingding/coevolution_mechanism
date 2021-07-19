# classifies each read into different categories for a fasta files
# see macpAndClassify300read() in mapClassPe

from mutTools import rev_complement, fasta_iter
from configx39 import expToInd, expToSeq
import mapClassPe as mcp
from constants import AT_SINGLE_FILTERED_DIR_M1, AT_SINGLE_FILTERED_DIR_M2, AT_SINGLE_CALLED_DIR_M1, AT_SINGLE_CALLED_DIR_M2

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

fastaInDir1 = AT_SINGLE_FILTERED_DIR_M1
dout1 = AT_SINGLE_CALLED_DIR_M1
classifySamples(fastaInDir1, dout1)

fastaInDir2 =AT_SINGLE_FILTERED_DIR_M2
dout2 = AT_SINGLE_CALLED_DIR_M2
classifySamples(fastaInDir2, dout2)
