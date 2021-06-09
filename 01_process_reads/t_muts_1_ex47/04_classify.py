import mapClassPe as mcp
from mutTools import fasta_iter_py3
import yaml
import os
from os import listdir
from os.path import isfile, join
import pandas as pd

from pipelineTools import class_samples

#read config files
config_dics = yaml.safe_load(open('./ex47_config.yaml', 'r'))
df_config = pd.read_csv(open('./ex47_config.csv', 'r'))

if __name__ == '__main__':
    #for hiseq run 1
    #din = '/n/groups/marks/users/david/ex47/04_split_by_ats/'
    #dout = '/n/groups/marks/users/david/ex47/05_class/'


    #for hiseq run 2
    din = '/n/groups/marks/users/david/ex47/04_split_by_ats_3786W/'
    dout = '/n/groups/marks/users/david/ex47/05_class_3786W/'

    class_samples(din, dout, df_config)
