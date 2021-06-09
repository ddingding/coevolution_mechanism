# this script merged paired end reads using usearch
# and does quality filtering of the merged reads, see below.

# this is applicable to 2 MiSeq .bam files:
#machine 1 files are:



import os
from os import listdir
from os.path import isfile, join

from configx39 import expToInd
from mutTools import rev_complement


#for each pair of paired end reads (separate files)

def merge_reads(fastqPath, outputDir):
	with open(outputDir+'cmdsUsed.txt', 'w') as fout:
		for ind in expToInd.values():
			ind = rev_complement(ind)
			usearchCmd ='/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex8/Illumina/fastqProcessing/usearch '\
			+ '-fastq_mergepairs %s  -reverse %s -fastqout %s_merged.fastq' \
			% (fastqPath+ind+'_1_.fastq', fastqPath+ind+'_2_.fastq',outputDir+ind)
			print usearchCmd
			os.system(usearchCmd)
			fout.write(usearchCmd)

#merge from first miseq
fastqPath = '/Users/davidd/non_dropbox/3272/'
outputDir = '/Users/davidd/non_dropbox/3272_merged/'
merge_reads(fastqPath, outputDir)

#merge from second miseq
fastqPath = '/Users/davidd/non_dropbox/3252/'
outputDir = '/Users/davidd/non_dropbox/3252_merged/'
merge_reads(fastqPath, outputDir)

#2 using usearch to filter and truncate reads by quality score

#for the first set

def filter_quality(fastqPath, outputDir):
	fastq = [f for f in listdir(fastqPath) if isfile(join(fastqPath, f))]
	with open(outputDir+'cmdsUsed.txt', 'w') as fout:
		for fi in fastq:
			usearchCmd ='/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex8/Illumina/fastqProcessing/usearch '\
			+ '-fastq_filter %s -fastq_truncqual 20 -fastq_maxns 3 -fastq_maxee 0.5 -fastq_ascii 33 -fastaout %s.fasta' \
			% (fastqPath+fi, outputDir+fi[:-6])

			print usearchCmd
			os.system(usearchCmd)
			fout.write(usearchCmd)

fastqPath1= '/Users/davidd/non_dropbox/3252_merged/'
outputDir1 = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/03fastqFilteredMi1/'
filter_quality(fastqPath1, outputDir1)

fastqPath2= '/Users/davidd/non_dropbox/3272_merged/'
outputDir2 ='/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/03fastqFilteredMi2/'
filter_quality(fastqPath2, outputDir2)
