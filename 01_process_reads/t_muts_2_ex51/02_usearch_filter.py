# 2 using usearch to filter and truncate reads by quality score

from os import listdir
from os.path import isfile, join
from pipelineTools import filter_fastq_quality
from constants import T_SINGLE_2_MERGED_DIR_M1, T_SINGLE_2_MERGED_DIR_M2, T_SINGLE_2_FILTERED_DIR_M1, \
    T_SINGLE_2_FILTERED_DIR_M2, VSEARCH_PATH


def filter_all_fastq(fastqPath, outputDir):
    fastq = [
        f
        for f in listdir(fastqPath)
        if isfile(join(fastqPath, f)) and f.endswith("extendedFrags.fastq")
    ]

    print("filtering files ", fastq)
    with open(outputDir + "cmdsUsed.txt", "w") as fout:
        for fi in fastq:
            filter_cmd = filter_fastq_quality(fastqPath + fi, outputDir + fi[:-6], vsearch_path = VSEARCH_PATH)
            fout.write(filter_cmd + "\n")


# for hiseq run 1
fastqPath = T_SINGLE_2_MERGED_DIR_M1
outputDir = T_SINGLE_2_FILTERED_DIR_M1
filter_all_fastq(fastqPath, outputDir)

# for hiseq run 2
fastqPath = T_SINGLE_2_MERGED_DIR_M2
outputDir = T_SINGLE_2_FILTERED_DIR_M2
filter_all_fastq(fastqPath, outputDir)
