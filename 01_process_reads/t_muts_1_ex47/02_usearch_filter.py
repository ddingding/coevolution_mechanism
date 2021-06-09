import os
from os import listdir
from os.path import isfile, join
from pipelineTools import filter_fastq_quality

#2 using usearch to filter and truncate reads by quality score



if __name__ == '__main__':
	#find all files in one dir and truncate and write to outputDir

	#for hiseq run 1
	fastqPath = '/n/groups/marks/users/david/ex47/02merged/'
	outputDir = '/n/groups/marks/users/david/ex47/03filtered/'

	#for hiseq run 2
	#fastqPath = '/n/groups/marks/users/david/ex47/02merged_2_3786W/'
	#outputDir = '/n/groups/marks/users/david/ex47/03filtered_2_3786W/'

	done_fastq_nums = [f.rstrip('_merged_fasta')[-3:] for f in listdir(outputDir) if
				isfile(join(outputDir, f)) and f.endswith('merged.fasta')]

	fastq = [f for f in listdir(fastqPath) if isfile(join(fastqPath, f)) and
			f.endswith('merged.fastq') and
			f.rstrip('-_merged.fastq')[-3:] not in done_fastq_nums]

	print('filtering files ', fastq)
	with open(outputDir+'cmdsUsed.txt', 'w') as fout:
		for fi in fastq:
			filter_cmd = filter_fastq_quality(fastqPath+fi, outputDir+fi[:-6])
			fout.write(filter_cmd +'\n')
