# classifies each read into different categories for a fasta files
# see macpAndClassify300read() in mapClassPe

from mutTools import fasta_iter_py3
import mapClassPe as mcp
from constants import (
    AT_SINGLE_FILTERED_DIR_M1,
    AT_SINGLE_FILTERED_DIR_M2,
    AT_SINGLE_CALLED_DIR_M1,
    AT_SINGLE_CALLED_DIR_M2,
)
from os import listdir


def classifySamples(fastaInDir, dout, template="pard"):
    print("classifying..." + fastaInDir)
    fastas = [f for f in listdir(fastaInDir) if f.endswith(".fasta")]
    for f in fastas:
        with open(dout + f[: -len("_merged.fasta")] + "_class.csv", "w") as fout:
            c = 0
            for rec in fasta_iter_py3(fastaInDir + f):
                writeL = mcp.mapAndClassify300Read(rec.header, rec.sequence, template)
                fout.write("\t".join(writeL) + "\n")
                c += 1
                if c % 10000 == 1:
                    print(c)


# classifying machine 1 reads
fastaInDir1 = AT_SINGLE_FILTERED_DIR_M1
dout1 = AT_SINGLE_CALLED_DIR_M1
classifySamples(fastaInDir1, dout1)

# classifying machine 2 reads
fastaInDir2 = AT_SINGLE_FILTERED_DIR_M2
dout2 = AT_SINGLE_CALLED_DIR_M2
classifySamples(fastaInDir2, dout2)
