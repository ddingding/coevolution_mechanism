# split reads by their antitoxin barcode (4 nucleotide at C-terminus of toxin gene)

from os import listdir
from os.path import isfile, join
import yaml
import pandas as pd
from pipelineTools import demultiplex_fastas, map_bmc_to_primer

from constants import T_SINGLE_2_FILTERED_DIR_M1, T_SINGLE_2_FILTERED_DIR_M2, T_SINGLE_2_AT_SPLIT_DIR_M1, T_SINGLE_2_AT_SPLIT_DIR_M2
################################################################################
# read config files
config_dics = yaml.safe_load(open("./ex51_config.yaml", "r"))
df_config = pd.read_csv(open("./ex51_config.csv", "r"))
################################################################################
# specify directories of files to process

def split_all_files(fasta_dir_in, fasta_dout):
    # select fasta files to process
    fastas = [
        fasta_dir_in + "190716Lau_" + f.split("_")[1]
        for f in listdir(fasta_dir_in)
        if isfile(join(fasta_dir_in, f)) and f.endswith(".fasta")
    ]
    '''
    # all the files with a _split_stats.csv are done files.
    done_prim_nums = [
        f[:3]
        for f in listdir(fasta_dout)
        if isfile(join(fasta_dout, f))
        if f.endswith("_split_stats.csv")
    ]
    '''
    print(fastas)
    # take only fasta files which are not done yet.

    fastas = [
        f
        for f in fastas
        #if map_bmc_to_primer(config_dics, f.split("_")[-1][:-2]) not in done_prim_nums
    ]

    print(fastas)

    demultiplex_fastas(fastas, fasta_dout, df_config, config_dics, exp_num=2, execute=True)

# hiseq 1 run
fasta_dir_in = T_SINGLE_2_FILTERED_DIR_M1
fasta_dout = T_SINGLE_2_AT_SPLIT_DIR_M1
split_all_files(fasta_dir_in, fasta_dout)

# hiseq 2 run
fasta_dir_in = T_SINGLE_2_FILTERED_DIR_M2
fasta_dout = T_SINGLE_2_AT_SPLIT_DIR_M2
split_all_files(fasta_dir_in, fasta_dout)

