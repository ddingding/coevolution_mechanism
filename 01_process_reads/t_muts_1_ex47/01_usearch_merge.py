## 1() filter by phred score.
# using edgar's fastq filter usearch
# switched to vsearch to handle larger files ~ >5gb

import os
from os import listdir
from os.path import isfile, join
from pipelineTools import merge_paired_reads_vsearch_o2

################################################################################

################################################################################

if __name__ == '__main__':
	#merge dirs

	#for the first machine:
	fastqPath = '/n/groups/marks/users/david/ex47/01fastq/'
	outputDir = '/n/groups/marks/users/david/ex47/02merged/'

	# for the second machine
	#fastqPath = '/n/groups/marks/users/david/ex47/01fastq_2_3786W/'
	#outputDir = '/n/groups/marks/users/david/ex47/02merged_2_3786W/'


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
