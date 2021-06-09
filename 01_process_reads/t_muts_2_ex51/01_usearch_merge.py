## merging paired end reads together.

import os
from os import listdir
from os.path import isfile, join
from pipelineTools import merge_paired_reads_vsearch_o2

################################################################################

################################################################################

if __name__ == "__main__":
    # merge dirs

    # for the first machine:
    fastqPath = "/n/groups/marks/users/david/ex51/00fastq_raw/"
    outputDir = "/n/groups/marks/users/david/ex51/02merged/"

    # for the second machine
    # fastqPath = '/n/groups/marks/users/david/ex51/00fastq_raw_4021W/'
    # outputDir = '/n/groups/marks/users/david/ex51/02merged_4021W/'

    n_to_do = 1000
    c = 0

    with open(outputDir + "cmdsUsed.txt", "w") as fout:
        # get a list of files that are already merged.
        done_nums = [
            f.rstrip("-_merged.fastq")[-3:]
            for f in listdir(outputDir)
            if isfile(join(outputDir, f))
        ]
        print("doing filepahts:{}{}".format(fastqPath, outputDir))
        print("file numbers already done:", done_nums)
        fastq = [
            f.rstrip("_sequence.fastq")[:-1]
            for f in listdir(fastqPath)
            if isfile(join(fastqPath, f))
        ]  # and
        # f.rstrip('_sequence.fastq')[-5:-2] not in done_nums]
        # fastq = ['190716Lau_D19-7910-2_']
        print("fastq files to do ", sorted(fastq))

        # fs

        for fa in fastq:
            if c < n_to_do:

                flash_cmd = "/n/groups/marks/users/david/apps/FLASH-1.2.11/flash {}{}1_sequence.fastq {}{}2_sequence.fastq -o {} -d {}".format(
                    fastqPath, fa, fastqPath, fa, fa, outputDir
                )
                print(flash_cmd)
                os.system(flash_cmd)
                # print('merging files: {} \n {}'.format(fastqPath+fa+'1_sequence.fastq', fastqPath+fa+'2_sequence.fastq'))
                c += 1
