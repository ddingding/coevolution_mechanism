# 2 using vsearch to filter and truncate reads by quality score
from os import listdir
from os.path import isfile, join
from pipelineTools import filter_fastq_quality
from constants import T_SINGLE_1_MERGED_DIR_M1, T_SINGLE_1_MERGED_DIR_M2, T_SINGLE_1_FILTERED_DIR_M1, \
    T_SINGLE_1_FILTERED_DIR_M2


def filter_fastq_quality_all(fastqPath, outputDir):
    # fetch all files to quality-filter, and filter them into a new output directory using vsearch
    done_fastq_nums = [f.rstrip('_merged_fasta')[-3:] for f in listdir(outputDir) if
                       isfile(join(outputDir, f)) and f.endswith('merged.fasta')]

    fastq = [f for f in listdir(fastqPath) if isfile(join(fastqPath, f)) and
             f.endswith('merged.fastq') and
             f.rstrip('-_merged.fastq')[-3:] not in done_fastq_nums]

    print('filtering files ', fastq)
    with open(outputDir + 'cmdsUsed.txt', 'w') as fout:
        for fi in fastq:
            filter_cmd = filter_fastq_quality(fastqPath + fi, outputDir + fi[:-6])
            fout.write(filter_cmd + '\n')


# for hiseq run 1
fastqPath = T_SINGLE_1_MERGED_DIR_M1
outputDir = T_SINGLE_1_FILTERED_DIR_M1
filter_fastq_quality_all(fastqPath, outputDir)

# for hiseq run 2
fastqPath = T_SINGLE_1_MERGED_DIR_M2
outputDir = T_SINGLE_1_FILTERED_DIR_M2
filter_fastq_quality_all(fastqPath, outputDir)
