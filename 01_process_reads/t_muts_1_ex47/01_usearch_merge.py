# call vsearch to merge paired end reads.

from os import listdir
from os.path import isfile, join
from pipelineTools import merge_paired_reads_vsearch_o2
from constants import T_SINGLE_1_CONCAT_DIR_M1, T_SINGLE_1_CONCAT_DIR_M2, T_SINGLE_1_MERGED_DIR_M1, T_SINGLE_1_MERGED_DIR_M2


def merge_reads(fastqPath, outputDir):
	with open(outputDir+'cmdsUsed.txt', 'w') as fout:
		# get a list of files that are already merged.
		done_nums = [f.rstrip('-_merged.fastq')[-3:] for f in listdir(outputDir)
					if isfile(join(outputDir, f))]
		print('file numbers already done:',done_nums)
		fastq = [f.rstrip('_sequence.fastq')[:-1]
				for f in listdir(fastqPath) if isfile(join(fastqPath, f)) and
				f.rstrip('_sequence.fastq')[-5:-2] not in done_nums]
		print('fastq files to do ',sorted(fastq))

		#fs should end in 190311Lau_D19-2172-
		for fa in fastq:
			vsearch_cmd = merge_paired_reads_vsearch_o2(fastqPath+fa+'1_sequence.fastq',
								fastqPath+fa+'2_sequence.fastq',outputDir+fa[:-1])

			fout.write(vsearch_cmd + '\n')

#for the first machine:
fastqPath = T_SINGLE_1_CONCAT_DIR_M1
outputDir = T_SINGLE_1_MERGED_DIR_M1
merge_reads(fastqPath, outputDir)

# for the second machine
fastqPath = T_SINGLE_1_CONCAT_DIR_M2
outputDir = T_SINGLE_1_MERGED_DIR_M2
merge_reads(fastqPath, outputDir)

