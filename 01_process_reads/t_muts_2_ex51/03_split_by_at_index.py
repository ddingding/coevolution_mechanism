# split reads by their antitoxin barcode (4 nucleotide at C-terminus of toxin gene)

from os import listdir
from os.path import isfile, join
import yaml
import pandas as pd
from pipelineTools import demultiplex_fastas, map_bmc_to_primer

################################################################################
# read config files
config_dics = yaml.safe_load(open("./ex51_config.yaml", "r"))
df_config = pd.read_csv(open("./ex51_config.csv", "r"))
################################################################################
# specify directories of files to process

def split_all_files(fasta_dir_in, fasta_dout):
    # select fasta files to process
    fastas = [
        fasta_dir_in + "190716Lau_" + f.split("_")[-2]
        for f in listdir(fasta_dir_in)
        if isfile(join(fasta_dir_in, f)) and f.endswith(".fasta")
    ]

    # all the files with a _split_stats.csv are done files.
    done_prim_nums = [
        f[:3]
        for f in listdir(fasta_dout)
        if isfile(join(fasta_dout, f))
        if f.endswith("_split_stats.csv")
    ]

    print(fastas)
    # take only fasta files which are not done yet.

    # changed this for new filenames
    fastas = [
        f
        for f in fastas
        if map_bmc_to_primer(config_dics, f.split("_")[-1][:-2]) not in done_prim_nums
    ]

    print(fastas)

    demultiplex_fastas(fastas, fasta_dout, df_config, config_dics)

# hiseq 1 run
fasta_dir_in = "/n/groups/marks/users/david/ex51/03filtered/"
fasta_dout = "/n/groups/marks/users/david/ex51/04_split_by_ats/"
split_all_files(fasta_dir_in, fasta_dout)

# hiseq 2 run
fasta_dir_in = "/n/groups/marks/users/david/ex51/03filtered_4021W/"
fasta_dout = "/n/groups/marks/users/david/ex51/04_split_by_ats_4021W/"
split_all_files(fasta_dir_in, fasta_dout)

