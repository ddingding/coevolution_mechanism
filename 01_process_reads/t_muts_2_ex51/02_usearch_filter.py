# 2 using usearch to filter and truncate reads by quality score

from os import listdir
from os.path import isfile, join
from pipelineTools import filter_fastq_quality

def filter_all_fastq(fastqPath, outputDir):
    fastq = [
        f
        for f in listdir(fastqPath)
        if isfile(join(fastqPath, f)) and f.endswith("extendedFrags.fastq")
    ]

    print("filtering files ", fastq)
    with open(outputDir + "cmdsUsed.txt", "w") as fout:
        for fi in fastq:
            filter_cmd = filter_fastq_quality(fastqPath + fi, outputDir + fi[:-6])
            fout.write(filter_cmd + "\n")

# for hiseq run 1
fastqPath = "/n/groups/marks/users/david/ex51/02merged/"
outputDir = "/n/groups/marks/users/david/ex51/03filtered/"
filter_all_fastq(fastqPath, outputDir)

# for hiseq run 2
fastqPath = '/n/groups/marks/users/david/ex51/02merged_4021W/'
outputDir = '/n/groups/marks/users/david/ex51/03filtered_4021W/'
filter_all_fastq(fastqPath, outputDir)