'''
## using flash 1.2.11 to merge paired end reads
the sample identified is in the middle of the filename: 'D21-169', refer to ex58_config.csv for description of that
sample
each paired end read are split in the files ending in '1_sequence.fastq' and '2_sequence.fastq'
like so:
210105Lau_D21-169_1_sequence.fastq	210105Lau_D21-169_2_sequence.fastq
210105Lau_D21-170_1_sequence.fastq	210105Lau_D21-170_2_sequence.fastq
210105Lau_D21-171_1_sequence.fastq	210105Lau_D21-171_2_sequence.fastq
210105Lau_D21-172_1_sequence.fastq	210105Lau_D21-172_2_sequence.fastq
210105Lau_D21-173_1_sequence.fastq	210105Lau_D21-173_2_sequence.fastq
210105Lau_D21-174_1_sequence.fastq	210105Lau_D21-174_2_sequence.fastq
210105Lau_D21-175_1_sequence.fastq	210105Lau_D21-175_2_sequence.fastq
210105Lau_D21-176_1_sequence.fastq	210105Lau_D21-176_2_sequence.fastq
210105Lau_D21-177_1_sequence.fastq	210105Lau_D21-177_2_sequence.fastq
210105Lau_D21-178_1_sequence.fastq	210105Lau_D21-178_2_sequence.fastq
210105Lau_D21-179_1_sequence.fastq	210105Lau_D21-179_2_sequence.fastq
210105Lau_D21-180_1_sequence.fastq	210105Lau_D21-180_2_sequence.fastq
210105Lau_D21-181_1_sequence.fastq	210105Lau_D21-181_2_sequence.fastq
210105Lau_D21-182_1_sequence.fastq	210105Lau_D21-182_2_sequence.fastq
210105Lau_D21-183_1_sequence.fastq	210105Lau_D21-183_2_sequence.fastq
210105Lau_D21-184_1_sequence.fastq	210105Lau_D21-184_2_sequence.fastq
210105Lau_D21-185_1_sequence.fastq	210105Lau_D21-185_2_sequence.fastq
210105Lau_D21-186_1_sequence.fastq	210105Lau_D21-186_2_sequence.fastq
210105Lau_D21-187_1_sequence.fastq	210105Lau_D21-187_2_sequence.fastq
210105Lau_D21-188_1_sequence.fastq	210105Lau_D21-188_2_sequence.fastq
210105Lau_D21-189_1_sequence.fastq	210105Lau_D21-189_2_sequence.fastq
210105Lau_D21-190_1_sequence.fastq	210105Lau_D21-190_2_sequence.fastq
210105Lau_D21-191_1_sequence.fastq	210105Lau_D21-191_2_sequence.fastq
210105Lau_D21-192_1_sequence.fastq	210105Lau_D21-192_2_sequence.fastq
210105Lau_D21-193_1_sequence.fastq	210105Lau_D21-193_2_sequence.fastq
210105Lau_D21-194_1_sequence.fastq	210105Lau_D21-194_2_sequence.fastq
210105Lau_D21-195_1_sequence.fastq	210105Lau_D21-195_2_sequence.fastq
210105Lau_D21-196_1_sequence.fastq	210105Lau_D21-196_2_sequence.fastq
210105Lau_D21-197_1_sequence.fastq	210105Lau_D21-197_2_sequence.fastq
210105Lau_D21-198_1_sequence.fastq	210105Lau_D21-198_2_sequence.fastq
210105Lau_D21-199_1_sequence.fastq	210105Lau_D21-199_2_sequence.fastq
210105Lau_D21-200_1_sequence.fastq	210105Lau_D21-200_2_sequence.fastq
210105Lau_D21-201_1_sequence.fastq	210105Lau_D21-201_2_sequence.fastq
210105Lau_D21-202_1_sequence.fastq	210105Lau_D21-202_2_sequence.fastq
210105Lau_D21-203_1_sequence.fastq	210105Lau_D21-203_2_sequence.fastq
210105Lau_D21-204_1_sequence.fastq	210105Lau_D21-204_2_sequence.fastq
210105Lau_D21-205_1_sequence.fastq	210105Lau_D21-205_2_sequence.fastq
210105Lau_D21-206_1_sequence.fastq	210105Lau_D21-206_2_sequence.fastq
210105Lau_D21-207_1_sequence.fastq	210105Lau_D21-207_2_sequence.fastq
210105Lau_D21-208_1_sequence.fastq	210105Lau_D21-208_2_sequence.fastq
210105Lau_D21-209_1_sequence.fastq	210105Lau_D21-209_2_sequence.fastq
210105Lau_D21-210_1_sequence.fastq	210105Lau_D21-210_2_sequence.fastq
210105Lau_D21-211_1_sequence.fastq	210105Lau_D21-211_2_sequence.fastq
210105Lau_D21-212_1_sequence.fastq	210105Lau_D21-212_2_sequence.fastq
210105Lau_D21-213_1_sequence.fastq	210105Lau_D21-213_2_sequence.fastq
210105Lau_D21-214_1_sequence.fastq	210105Lau_D21-214_2_sequence.fastq
210105Lau_D21-215_1_sequence.fastq	210105Lau_D21-215_2_sequence.fastq
210105Lau_D21-216_1_sequence.fastq	210105Lau_D21-216_2_sequence.fastq
'''

import os
from os import listdir
from os.path import isfile, join
from constants import AT_COMBO_RAW_DIR, AT_COMBO_MERGED_DIR
fastqPath = AT_COMBO_RAW_DIR
outputDir = AT_COMBO_MERGED_DIR

# merge dirs
n_to_do = 1000
c = 0

with open(outputDir + 'cmdsUsed.txt', 'w') as fout:
    # get a list of files that are already merged.
    done_nums = [f.rstrip('-_merged.fastq')[-3:] for f in listdir(outputDir)
                 if isfile(join(outputDir, f))]
    print('doing filepahts:{}{}'.format(fastqPath, outputDir))
    print('file numbers already done:', done_nums)
    fastq = [f.rstrip('_sequence.fastq')[:-1]
             for f in listdir(fastqPath) if isfile(join(fastqPath, f))]
    print('fastq files to do ', sorted(fastq))

    # merge each file
    for fa in fastq:
        if c < n_to_do:
            flash_cmd = '/n/groups/marks/users/david/apps/FLASH-1.2.11/flash {}{}1_sequence.fastq {}{' \
						'}2_sequence.fastq -o {} -d {}'.format(
                fastqPath, fa, fastqPath, fa, fa, outputDir)
            print(flash_cmd)
            os.system(flash_cmd)
            # print('merging files: {} \n {}'.format(fastqPath+fa+'1_sequence.fastq', fastqPath+fa+'2_sequence.fastq'))
            c += 1
