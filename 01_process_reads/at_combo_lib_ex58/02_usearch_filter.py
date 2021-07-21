# using usearch to filter and truncate reads by quality score


from os import listdir
from os.path import isfile, join
from pipelineTools import filter_fastq_quality
from constants import AT_COMBO_MERGED_DIR, AT_COMBO_FILTERED_DIR, VSEARCH_PATH

fastqPath = AT_COMBO_MERGED_DIR
outputDir = AT_COMBO_FILTERED_DIR

fastq = [f for f in listdir(fastqPath) if isfile(join(fastqPath, f)) and
					f.endswith('extendedFrags.fastq')]

print('filtering files ', fastq)
with open(outputDir+'cmdsUsed.txt', 'w') as fout:
	for fi in fastq:
		filter_cmd = filter_fastq_quality(fastqPath+fi, outputDir+fi[:-6], vsearch_path = VSEARCH_PATH)
		fout.write(filter_cmd +'\n')
