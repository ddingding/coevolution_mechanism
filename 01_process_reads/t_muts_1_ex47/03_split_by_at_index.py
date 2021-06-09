import os
from os import listdir
from os.path import isfile, join
import yaml
from mutTools import fasta_iter_py3
import pandas as pd #pandas can be imported in the python3 environment on o2
from pipelineTools import demultiplex_fastas, map_bmc_to_primer

################################################################################
#read config files
config_dics = yaml.safe_load(open('./ex47_config.yaml', 'r'))
df_config = pd.read_csv(open('./ex47_config.csv', 'r'))
################################################################################
# specify directories of files to process

#hiseq 1 run
fasta_dir_in = '/n/groups/marks/users/david/ex47/03filtered/'
fasta_dout = '/n/groups/marks/users/david/ex47/04_split_by_ats/'

#hiseq 2
#fasta_dir_in = '/n/groups/marks/users/david/ex47/03filtered_2_3786W/'
#fasta_dout = '/n/groups/marks/users/david/ex47/04_split_by_ats_3786W/'
################################################################################

#############  MAIN ###########################################################

if __name__ == '__main__':
    #select fasta files to process

    fastas = [fasta_dir_in+f for f in listdir(fasta_dir_in) if isfile(join(fasta_dir_in, f)) and
                f.endswith('.fasta')]

    # all the files with a _split_stats.csv are done files.
    done_prim_nums = [f[:3] for f in listdir(fasta_dout) if isfile(join(fasta_dout, f))
                    if f.endswith('_split_stats.csv')]


    #take only fasta files which are not done yet.
    fastas = [f for f in fastas
        if map_bmc_to_primer(config_dics,f.split('_')[-2]) not in done_prim_nums]

    print(fastas)

    demultiplex_fastas(fastas, fasta_dout,df_config, config_dics)
