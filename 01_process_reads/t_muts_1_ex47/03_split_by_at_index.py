'''
Some of these illumina read count files contain toxin single variant library in the background of different antitoxin
mutant backgrounds.
This script deconvolutes toxin single mutants in different antitoxin variant backgrounds from the same file.
To read out the antitoxin single mutant backgrounds, there is a 4 nucleotide barcode at the 3' end of the toxin
vector (that was matched with a particular antitoxin background during transformation fo the toxin library into cells
containing that particular antitoxin variant)
see the ex47_config.yaml file for pairing of the barcodes to the antitoxin variant

'''

from os import listdir
from os.path import isfile, join
import yaml
import pandas as pd
from pipelineTools import demultiplex_fastas, map_bmc_to_primer
from constants import T_SINGLE_1_FILTERED_DIR_M1, T_SINGLE_1_FILTERED_DIR_M2, T_SINGLE_1_AT_SPLIT_DIR_M1, T_SINGLE_1_AT_SPLIT_DIR_M2
################################################################################
# read config files
config_dics = yaml.safe_load(open('./ex47_config.yaml', 'r'))
df_config = pd.read_csv(open('./ex47_config.csv', 'r'))
################################################################################
def demultiplex_ats_all_fastas(fasta_dir_in, fasta_dout, config_dics, df_config):
    # select fasta files to process
    fastas = [fasta_dir_in + f for f in listdir(fasta_dir_in) if isfile(join(fasta_dir_in, f)) and
              f.endswith('.fasta')]

    # all the files with a _split_stats.csv are done files.
    done_prim_nums = [f[:3] for f in listdir(fasta_dout) if isfile(join(fasta_dout, f))
                      if f.endswith('_split_stats.csv')]

    # take only fasta files which are not done yet.
    fastas = [f for f in fastas
              if map_bmc_to_primer(config_dics, f.split('_')[-2]) not in done_prim_nums]

    print(fastas)

    demultiplex_fastas(fastas, fasta_dout, df_config, config_dics)


# hiseq machine 1 results
fasta_dir_in = T_SINGLE_1_FILTERED_DIR_M1
fasta_dout = T_SINGLE_1_AT_SPLIT_DIR_M1
demultiplex_ats_all_fastas(fasta_dir_in, fasta_dout, config_dics, df_config)

# hiseq machine 2 results
fasta_dir_in = T_SINGLE_1_FILTERED_DIR_M2
fasta_dout = T_SINGLE_1_AT_SPLIT_DIR_M2
demultiplex_ats_all_fastas(fasta_dir_in, fasta_dout, config_dics, df_config)

