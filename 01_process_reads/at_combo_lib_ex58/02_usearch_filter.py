import os
from os import listdir
from os.path import isfile, join
from pipelineTools import filter_fastq_quality

#2 using usearch to filter and truncate reads by quality score



if __name__ == '__main__':
	#find all files in one dir and truncate and write to outputDir

	fastqPath = '/n/groups/marks/users/david/ex58/01_merged/'
	outputDir = '/n/groups/marks/users/david/ex58/02_filtered/'

	done_fastq_nums = [f.rstrip('_merged_fasta')[-3:] for f in listdir(outputDir) if
				isfile(join(outputDir, f)) and f.endswith('merged.fasta')]


	fastq = [f for f in listdir(fastqPath) if isfile(join(fastqPath, f)) and
                      	f.endswith('extendedFrags.fastq')]

	print('filtering files ', fastq)
	with open(outputDir+'cmdsUsed.txt', 'w') as fout:
		for fi in fastq:
			filter_cmd = filter_fastq_quality(fastqPath+fi, outputDir+fi[:-6])
			fout.write(filter_cmd +'\n')
