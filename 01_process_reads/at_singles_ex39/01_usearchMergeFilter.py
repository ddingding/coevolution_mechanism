'''
# this script merged paired end reads using usearch
# and does quality filtering of the merged reads, see below.

# this is applicable to 2 MiSeq files
miseq run 1 files are:
180716Lau_D18-6083_phiX_bestmap.fastq.gz
180716Lau_D18-6084_phiX_bestmap.fastq.gz
180716Lau_D18-6093_phiX_bestmap.fastq.gz
180716Lau_D18-6094_phiX_bestmap.fastq.gz

miseq run 2 files are:
180716Lau_D18-6083_phiX_bestmap_m2.fastq.gz
180716Lau_D18-6084_phiX_bestmap_m2.fastq.gz
180716Lau_D18-6093_phiX_bestmap_m2.fastq.gz
180716Lau_D18-6094_phiX_bestmap_m2.fastq.gz

Note: The paired end reads are both in each of these files on the SRA, and will have to be split into 2 separate
files to be compatible with this script.
'''
import os
from os import listdir
from os.path import isfile, join

from configx39 import expToInd
from mutTools import rev_complement
from constants import AT_SINGLE_RAW_DIR_M1, AT_SINGLE_RAW_DIR_M2, AT_SINGLE_MERGED_DIR_M1, AT_SINGLE_MERGED_DIR_M2, AT_SINGLE_FILTERED_DIR_M1, AT_SINGLE_FILTERED_DIR_M2

########################################################################################################################
# merge each pair of paired end reads from separate files
def merge_reads(fastqPath, outputDir):
    with open(outputDir + 'cmdsUsed.txt', 'w') as fout:
        for ind in expToInd.values():
            ind = rev_complement(ind)
            usearchCmd = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex8/Illumina/fastqProcessing' \
						 '/usearch ' \
                         + '-fastq_mergepairs %s  -reverse %s -fastqout %s_merged.fastq' \
                         % (fastqPath + ind + '_1_.fastq', fastqPath + ind + '_2_.fastq', outputDir + ind)
            print usearchCmd
            os.system(usearchCmd)
            fout.write(usearchCmd)

# merge from first miseq
fastqPath = AT_SINGLE_RAW_DIR_M1
outputDir = AT_SINGLE_MERGED_DIR_M1
merge_reads(fastqPath, outputDir)

# merge from second miseq
fastqPath = AT_SINGLE_RAW_DIR_M2
outputDir = AT_SINGLE_MERGED_DIR_M2
merge_reads(fastqPath, outputDir)

########################################################################################################################
# 2 using usearch to filter and truncate reads by quality score

# for the first set

def filter_quality(fastqPath, outputDir):
    fastq = [f for f in listdir(fastqPath) if isfile(join(fastqPath, f))]
    with open(outputDir + 'cmdsUsed.txt', 'w') as fout:
        for fi in fastq:
            usearchCmd = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex8/Illumina/fastqProcessing' \
						 '/usearch ' \
                         + '-fastq_filter %s -fastq_truncqual 20 -fastq_maxns 3 -fastq_maxee 0.5 -fastq_ascii 33 ' \
						   '-fastaout %s.fasta' \
                         % (fastqPath + fi, outputDir + fi[:-6])

            print usearchCmd
            os.system(usearchCmd)
            fout.write(usearchCmd)


fastqPath1 = AT_SINGLE_MERGED_DIR_M1
outputDir1 = AT_SINGLE_FILTERED_DIR_M1
filter_quality(fastqPath1, outputDir1)

fastqPath2 = AT_SINGLE_MERGED_DIR_M2
outputDir2 = AT_SINGLE_FILTERED_DIR_M2
filter_quality(fastqPath2, outputDir2)
