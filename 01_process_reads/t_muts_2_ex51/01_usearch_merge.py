'''
## merging paired end reads together using flash.

Machine 1 paired end reads are (with corresponding paired end read files in a particular line):
190716Lau_D19-7910-1_1_sequence.fastq	190716Lau_D19-7910-1_2_sequence.fastq
190716Lau_D19-7911-2_1_sequence.fastq	190716Lau_D19-7911-2_2_sequence.fastq
190716Lau_D19-7912-2_1_sequence.fastq	190716Lau_D19-7912-2_2_sequence.fastq
190716Lau_D19-7913-1_1_sequence.fastq	190716Lau_D19-7913-1_2_sequence.fastq
190716Lau_D19-7914-1_1_sequence.fastq	190716Lau_D19-7914-1_2_sequence.fastq
190716Lau_D19-7915-2_1_sequence.fastq	190716Lau_D19-7915-2_2_sequence.fastq
190716Lau_D19-7916-2_1_sequence.fastq	190716Lau_D19-7916-2_2_sequence.fastq
190716Lau_D19-7917-2_1_sequence.fastq	190716Lau_D19-7917-2_2_sequence.fastq
190716Lau_D19-7918-1_1_sequence.fastq	190716Lau_D19-7918-1_2_sequence.fastq
190716Lau_D19-7919-1_1_sequence.fastq	190716Lau_D19-7919-1_2_sequence.fastq
190716Lau_D19-7920-1_1_sequence.fastq	190716Lau_D19-7920-1_2_sequence.fastq
190716Lau_D19-7921-2_1_sequence.fastq	190716Lau_D19-7921-2_2_sequence.fastq
190716Lau_D19-7922-2_1_sequence.fastq	190716Lau_D19-7922-2_2_sequence.fastq
190716Lau_D19-7923-2_1_sequence.fastq	190716Lau_D19-7923-2_2_sequence.fastq
190716Lau_D19-7924-1_1_sequence.fastq	190716Lau_D19-7924-1_2_sequence.fastq
190716Lau_D19-7925-2_1_sequence.fastq	190716Lau_D19-7925-2_2_sequence.fastq
190716Lau_D19-7926-2_1_sequence.fastq	190716Lau_D19-7926-2_2_sequence.fastq
190716Lau_D19-7927-2_1_sequence.fastq	190716Lau_D19-7927-2_2_sequence.fastq
190716Lau_D19-7928-1_1_sequence.fastq	190716Lau_D19-7928-1_2_sequence.fastq
190716Lau_D19-7929-2_1_sequence.fastq	190716Lau_D19-7929-2_2_sequence.fastq
190716Lau_D19-7930-2_1_sequence.fastq	190716Lau_D19-7930-2_2_sequence.fastq
190716Lau_D19-7931-1_1_sequence.fastq	190716Lau_D19-7931-1_2_sequence.fastq
190716Lau_D19-7932-2_1_sequence.fastq	190716Lau_D19-7932-2_2_sequence.fastq
190716Lau_D19-7933-1_1_sequence.fastq	190716Lau_D19-7933-1_2_sequence.fastq
190716Lau_D19-7934-1_1_sequence.fastq	190716Lau_D19-7934-1_2_sequence.fastq
190716Lau_D19-7935-1_1_sequence.fastq	190716Lau_D19-7935-1_2_sequence.fastq
190716Lau_D19-7936-2_1_sequence.fastq	190716Lau_D19-7936-2_2_sequence.fastq
190716Lau_D19-7937-2_1_sequence.fastq	190716Lau_D19-7937-2_2_sequence.fastq

Machine 2 paired end reads are:
190716Lau_D19-7910-1_1_sequence_m2.fastq	190716Lau_D19-7910-1_2_sequence_m2.fastq
190716Lau_D19-7911-2_1_sequence_m2.fastq	190716Lau_D19-7911-2_2_sequence_m2.fastq
190716Lau_D19-7912-2_1_sequence_m2.fastq	190716Lau_D19-7912-2_2_sequence_m2.fastq
190716Lau_D19-7913-1_1_sequence_m2.fastq	190716Lau_D19-7913-1_2_sequence_m2.fastq
190716Lau_D19-7914-1_1_sequence_m2.fastq	190716Lau_D19-7914-1_2_sequence_m2.fastq
190716Lau_D19-7915-2_1_sequence_m2.fastq	190716Lau_D19-7915-2_2_sequence_m2.fastq
190716Lau_D19-7916-2_1_sequence_m2.fastq	190716Lau_D19-7916-2_2_sequence_m2.fastq
190716Lau_D19-7917-2_1_sequence_m2.fastq	190716Lau_D19-7917-2_2_sequence_m2.fastq
190716Lau_D19-7918-1_1_sequence_m2.fastq	190716Lau_D19-7918-1_2_sequence_m2.fastq
190716Lau_D19-7919-1_1_sequence_m2.fastq	190716Lau_D19-7919-1_2_sequence_m2.fastq
190716Lau_D19-7920-1_1_sequence_m2.fastq	190716Lau_D19-7920-1_2_sequence_m2.fastq
190716Lau_D19-7921-2_1_sequence_m2.fastq	190716Lau_D19-7921-2_2_sequence_m2.fastq
190716Lau_D19-7922-2_1_sequence_m2.fastq	190716Lau_D19-7922-2_2_sequence_m2.fastq
190716Lau_D19-7923-2_1_sequence_m2.fastq	190716Lau_D19-7923-2_2_sequence_m2.fastq
190716Lau_D19-7924-1_1_sequence_m2.fastq	190716Lau_D19-7924-1_2_sequence_m2.fastq
190716Lau_D19-7925-2_1_sequence_m2.fastq	190716Lau_D19-7925-2_2_sequence_m2.fastq
190716Lau_D19-7926-2_1_sequence_m2.fastq	190716Lau_D19-7926-2_2_sequence_m2.fastq
190716Lau_D19-7927-2_1_sequence_m2.fastq	190716Lau_D19-7927-2_2_sequence_m2.fastq
190716Lau_D19-7928-1_1_sequence_m2.fastq	190716Lau_D19-7928-1_2_sequence_m2.fastq
190716Lau_D19-7929-2_1_sequence_m2.fastq	190716Lau_D19-7929-2_2_sequence_m2.fastq
190716Lau_D19-7930-1_1_sequence_m2.fastq	190716Lau_D19-7930-1_2_sequence_m2.fastq
190716Lau_D19-7931-1_1_sequence_m2.fastq	190716Lau_D19-7931-1_2_sequence_m2.fastq
190716Lau_D19-7932-2_1_sequence_m2.fastq	190716Lau_D19-7932-2_2_sequence_m2.fastq
190716Lau_D19-7933-1_1_sequence_m2.fastq	190716Lau_D19-7933-1_2_sequence_m2.fastq
190716Lau_D19-7934-1_1_sequence_m2.fastq	190716Lau_D19-7934-1_2_sequence_m2.fastq
190716Lau_D19-7935-1_1_sequence_m2.fastq	190716Lau_D19-7935-1_2_sequence_m2.fastq
190716Lau_D19-7936-2_1_sequence_m2.fastq	190716Lau_D19-7936-2_2_sequence_m2.fastq
190716Lau_D19-7937-2_1_sequence_m2.fastq	190716Lau_D19-7937-2_2_sequence_m2.fastq
'''

import os
from os import listdir
from os.path import isfile, join

from constants import FLASH_PATH, T_SINGLE_2_RAW_DIR_M1, T_SINGLE_2_RAW_DIR_M2, T_SINGLE_2_MERGED_DIR_M1, \
    T_SINGLE_2_MERGED_DIR_M2


def merge_paired_end_reads(fastqPath, outputDir, fastq_suffix="_sequence.fastq",
                           flash_path='/n/groups/marks/users/david/apps/FLASH-1.2.11/flash'):
    n_to_do = 1000
    c = 0

    # get a list of files that are already merged.
    done_nums = [
        f.rstrip("-_merged.fastq")[-3:]
        for f in listdir(outputDir)
        if isfile(join(outputDir, f))
    ]
    print("doing filepahts:{}{}".format(fastqPath, outputDir))
    print("file numbers already done:", done_nums)
    fastq = [
        f.rstrip(fastq_suffix)[:-1]
        for f in listdir(fastqPath)
        if isfile(join(fastqPath, f))
    ]
    print("fastq files to do ", sorted(fastq))

    for fa in fastq:
        if c < n_to_do:
            flash_cmd = " {}{}1{} {}{}2{} -o {} -d {}".format(flash_path,
                                                              fastqPath, fa, fastq_suffix, fastqPath, fa, fastq_suffix,
                                                              fa, outputDir
                                                              )
            print(flash_cmd)
            os.system(flash_cmd)
            c += 1


# for the first machine:
fastqPath = T_SINGLE_2_RAW_DIR_M1
outputDir = T_SINGLE_2_MERGED_DIR_M1
merge_paired_end_reads(fastqPath, outputDir, fastq_suffix="_sequence.fastq", flash_path=FLASH_PATH)

# for the second machine
fastqPath = T_SINGLE_2_RAW_DIR_M2
outputDir = T_SINGLE_2_MERGED_DIR_M2
merge_paired_end_reads(fastqPath, outputDir, fastq_suffix="_sequence_m2.fastq", flash_path=FLASH_PATH)
