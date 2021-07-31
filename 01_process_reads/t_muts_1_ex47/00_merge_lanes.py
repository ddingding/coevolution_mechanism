'''
# merge files from 2 HISeq lanes together
for each machine, create a directory, and put these files in:
directory 1: Hiseq machine 1, lane 1 paired end reads (paired end read split across these 2 files):
190311Lau_D19-2172-1_1_sequence.fastq	190311Lau_D19-2172-1_2_sequence.fastq
190311Lau_D19-2173-1_1_sequence.fastq	190311Lau_D19-2173-1_2_sequence.fastq
190311Lau_D19-2174-1_1_sequence.fastq	190311Lau_D19-2174-1_2_sequence.fastq
190311Lau_D19-2175-1_1_sequence.fastq	190311Lau_D19-2175-1_2_sequence.fastq
190311Lau_D19-2176-1_1_sequence.fastq	190311Lau_D19-2176-1_2_sequence.fastq
190311Lau_D19-2177-1_1_sequence.fastq	190311Lau_D19-2177-1_2_sequence.fastq
190311Lau_D19-2178-1_1_sequence.fastq	190311Lau_D19-2178-1_2_sequence.fastq
190311Lau_D19-2179-1_1_sequence.fastq	190311Lau_D19-2179-1_2_sequence.fastq
190311Lau_D19-2180-1_1_sequence.fastq	190311Lau_D19-2180-1_2_sequence.fastq
190311Lau_D19-2181-1_1_sequence.fastq	190311Lau_D19-2181-1_2_sequence.fastq
190311Lau_D19-2182-1_1_sequence.fastq	190311Lau_D19-2182-1_2_sequence.fastq
190311Lau_D19-2183-1_1_sequence.fastq	190311Lau_D19-2183-1_2_sequence.fastq
190311Lau_D19-2184-1_1_sequence.fastq	190311Lau_D19-2184-1_2_sequence.fastq
190311Lau_D19-2185-1_1_sequence.fastq	190311Lau_D19-2185-1_2_sequence.fastq
190311Lau_D19-2186-1_1_sequence.fastq	190311Lau_D19-2186-1_2_sequence.fastq
190311Lau_D19-2187-1_1_sequence.fastq	190311Lau_D19-2187-1_2_sequence.fastq

directory 1: Hiseq machine 1, lane 2 paired end reads:
190311Lau_D19-2172-2_1_sequence.fastq	190311Lau_D19-2172-2_2_sequence.fastq
190311Lau_D19-2173-2_1_sequence.fastq	190311Lau_D19-2173-2_2_sequence.fastq
190311Lau_D19-2174-2_1_sequence.fastq	190311Lau_D19-2174-2_2_sequence.fastq
190311Lau_D19-2175-2_1_sequence.fastq	190311Lau_D19-2175-2_2_sequence.fastq
190311Lau_D19-2176-2_1_sequence.fastq	190311Lau_D19-2176-2_2_sequence.fastq
190311Lau_D19-2177-2_1_sequence.fastq	190311Lau_D19-2177-2_2_sequence.fastq
190311Lau_D19-2178-2_1_sequence.fastq	190311Lau_D19-2178-2_2_sequence.fastq
190311Lau_D19-2179-2_1_sequence.fastq	190311Lau_D19-2179-2_2_sequence.fastq
190311Lau_D19-2180-2_1_sequence.fastq	190311Lau_D19-2180-2_2_sequence.fastq
190311Lau_D19-2181-2_1_sequence.fastq	190311Lau_D19-2181-2_2_sequence.fastq
190311Lau_D19-2182-2_1_sequence.fastq	190311Lau_D19-2182-2_2_sequence.fastq
190311Lau_D19-2183-2_1_sequence.fastq	190311Lau_D19-2183-2_2_sequence.fastq
190311Lau_D19-2184-2_1_sequence.fastq	190311Lau_D19-2184-2_2_sequence.fastq
190311Lau_D19-2185-2_1_sequence.fastq	190311Lau_D19-2185-2_2_sequence.fastq
190311Lau_D19-2186-2_1_sequence.fastq	190311Lau_D19-2186-2_2_sequence.fastq
190311Lau_D19-2187-2_1_sequence.fastq	190311Lau_D19-2187-2_2_sequence.fastq

directory 2: Hiseq machine 2, lane 1 paired end reads:
190311Lau_D19-2172-1_1_sequence_m2.fastq	190311Lau_D19-2172-1_2_sequence_m2.fastq
190311Lau_D19-2173-1_1_sequence_m2.fastq	190311Lau_D19-2173-1_2_sequence_m2.fastq
190311Lau_D19-2174-1_1_sequence_m2.fastq	190311Lau_D19-2174-1_2_sequence_m2.fastq
190311Lau_D19-2175-1_1_sequence_m2.fastq	190311Lau_D19-2175-1_2_sequence_m2.fastq
190311Lau_D19-2176-1_1_sequence_m2.fastq	190311Lau_D19-2176-1_2_sequence_m2.fastq
190311Lau_D19-2177-1_1_sequence_m2.fastq	190311Lau_D19-2177-1_2_sequence_m2.fastq
190311Lau_D19-2178-1_1_sequence_m2.fastq	190311Lau_D19-2178-1_2_sequence_m2.fastq
190311Lau_D19-2179-1_1_sequence_m2.fastq	190311Lau_D19-2179-1_2_sequence_m2.fastq
190311Lau_D19-2180-1_1_sequence_m2.fastq	190311Lau_D19-2180-1_2_sequence_m2.fastq
190311Lau_D19-2181-1_1_sequence_m2.fastq	190311Lau_D19-2181-1_2_sequence_m2.fastq
190311Lau_D19-2182-1_1_sequence_m2.fastq	190311Lau_D19-2182-1_2_sequence_m2.fastq
190311Lau_D19-2183-1_1_sequence_m2.fastq	190311Lau_D19-2183-1_2_sequence_m2.fastq
190311Lau_D19-2184-1_1_sequence_m2.fastq	190311Lau_D19-2184-1_2_sequence_m2.fastq
190311Lau_D19-2185-1_1_sequence_m2.fastq	190311Lau_D19-2185-1_2_sequence_m2.fastq
190311Lau_D19-2186-1_1_sequence_m2.fastq	190311Lau_D19-2186-1_2_sequence_m2.fastq
190311Lau_D19-2187-1_1_sequence_m2.fastq	190311Lau_D19-2187-1_2_sequence_m2.fastq

directory 2: Hiseq machine 2, lane 2 paired end reads:
190311Lau_D19-2172-2_1_sequence_m2.fastq	190311Lau_D19-2172-2_2_sequence_m2.fastq
190311Lau_D19-2173-2_1_sequence_m2.fastq	190311Lau_D19-2173-2_2_sequence_m2.fastq
190311Lau_D19-2174-2_1_sequence_m2.fastq	190311Lau_D19-2174-2_2_sequence_m2.fastq
190311Lau_D19-2175-2_1_sequence_m2.fastq	190311Lau_D19-2175-2_2_sequence_m2.fastq
190311Lau_D19-2176-2_1_sequence_m2.fastq	190311Lau_D19-2176-2_2_sequence_m2.fastq
190311Lau_D19-2177-2_1_sequence_m2.fastq	190311Lau_D19-2177-2_2_sequence_m2.fastq
190311Lau_D19-2178-2_1_sequence_m2.fastq	190311Lau_D19-2178-2_2_sequence_m2.fastq
190311Lau_D19-2179-2_1_sequence_m2.fastq	190311Lau_D19-2179-2_2_sequence_m2.fastq
190311Lau_D19-2180-2_1_sequence_m2.fastq	190311Lau_D19-2180-2_2_sequence_m2.fastq
190311Lau_D19-2181-2_1_sequence_m2.fastq	190311Lau_D19-2181-2_2_sequence_m2.fastq
190311Lau_D19-2182-2_1_sequence_m2.fastq	190311Lau_D19-2182-2_2_sequence_m2.fastq
190311Lau_D19-2183-2_1_sequence_m2.fastq	190311Lau_D19-2183-2_2_sequence_m2.fastq
190311Lau_D19-2184-2_1_sequence_m2.fastq	190311Lau_D19-2184-2_2_sequence_m2.fastq
190311Lau_D19-2185-2_1_sequence_m2.fastq	190311Lau_D19-2185-2_2_sequence_m2.fastq
190311Lau_D19-2186-2_1_sequence_m2.fastq	190311Lau_D19-2186-2_2_sequence_m2.fastq
190311Lau_D19-2187-2_1_sequence_m2.fastq	190311Lau_D19-2187-2_2_sequence_m2.fastq
'''

from os import listdir
from os.path import isfile, join

from pipelineTools import concat_files
from constants import T_SINGLE_1_RAW_DIR_M1, T_SINGLE_1_RAW_DIR_M2, T_SINGLE_1_CONCAT_DIR_M1, T_SINGLE_1_CONCAT_DIR_M2

def merge_fwd_and_rev_reads(fastq, fastqPath, mergedPath, f_suffix = 'sequence.fastq'):
    # merges paired-end read files together for each sample
    for f in fastq:
        print('concatenating reads from 2 lanes for ', f)
        # merge fwd read 1 from lane 1 and lane 2
        concat_files([fastqPath + f + '1_1_'+f_suffix,
                      fastqPath + f + '2_1_'+ f_suffix],
                     mergedPath + f + '1_sequence.fastq')

        # merge rev read 2 from lane 1 and lane 2
        concat_files([fastqPath + f + '1_2_'+f_suffix,
                      fastqPath + f + '2_2_'+ f_suffix],
                     mergedPath + f + '2_sequence.fastq')

#for the first HiSeq machine reads
fastqPath = T_SINGLE_1_RAW_DIR_M1
mergedPath = T_SINGLE_1_CONCAT_DIR_M1
# get a list of unique sample name prefixes, like '190311Lau_D19-2172-'
fastq = set([f[:19] for f in listdir(fastqPath)
            if isfile(join(fastqPath, f)) and f.endswith('.fastq')])
print(fastq)
merge_fwd_and_rev_reads(fastq, fastqPath, mergedPath, f_suffix = 'sequence.fastq')

#for the second hiseq machine reads
fastqPath = T_SINGLE_1_RAW_DIR_M2
mergedPath = T_SINGLE_1_CONCAT_DIR_M2

fastq = set([f[:19] for f in listdir(fastqPath)
            if isfile(join(fastqPath, f)) and f.endswith('.fastq')])
print(fastq)
merge_fwd_and_rev_reads(fastq, fastqPath, mergedPath, f_suffix = 'sequence_m2.fastq')

