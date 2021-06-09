# merge files from 2 HISeq lanes together

import os
from os import listdir
from os.path import isfile, join
from pipelineTools import merge_files

#for the first HiSeq machine reads
fastqPath = '/n/groups/marks/users/david/ex47/00fastq_raw/'
mergedPath = '/n/groups/marks/users/david/ex47/01fastq/'


#for the second hiseq machine reads
#fastqPath = '/n/groups/marks/users/david/ex47/00fastq_raw_machine2_3786W/'
#mergedPath = '/n/groups/marks/users/david/ex47/01fastq_2/'

################################################################################


################################################################################
# do the first indices
select_bmc_indices_to_do = ['172','173','180','181']

fastq = set([f.rstrip('_sequence.fastq')[:-3] for f in listdir(fastqPath)
            if isfile(join(fastqPath, f)) and f.endswith('sequence.fastq')])
#fs end in 190311Lau_D19-2177-

#only do the files
if select_bmc_indices_to_do:
    fastq = [f for f in fastq if f[-4:-1] not in select_bmc_indices_to_do]

#encoding should be:
#190311Lau_D19-2172-1_1_sequence.fastq
#the first number stands for the lane number
#and the second number 1 stands for the forward or reverse read.

print(fastq)

for f in fastq:
    print('concatenating reads from 2 lanes for ',f)
    #merge fwd read 1 from lane 1 and lane 2
    concat_files([fastqPath + f + '1_1_sequence.fastq',
                fastqPath + f + '2_1_sequence.fastq'],
                mergedPath+f+'1_sequence.fastq')

    #merge rev read 2 from lane 1 and lane 2
    cocnat_files([fastqPath + f + '1_2_sequence.fastq',
                fastqPath + f + '2_2_sequence.fastq'],
                mergedPath+f+'2_sequence.fastq')
