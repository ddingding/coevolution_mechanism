# classify reads and call single mutants for the no antitoxin samples


import os
from os import listdir
import yaml
import pandas as pd
from pipelineTools import class_samples
from constants import T_SINGLE_2_FILTERED_DIR_M1, T_SINGLE_2_FILTERED_DIR_NO_AT_M1, T_SINGLE_2_FILTERED_DIR_M2, T_SINGLE_2_FILTERED_DIR_NO_AT_M2, T_SINGLE_2_CALLED_DIR_NO_AT_M1, T_SINGLE_2_CALLED_DIR_NO_AT_M2


# read config files
config_dics = yaml.safe_load(open("./ex51_config.yaml", "r"))
df_config = pd.read_csv(open("./ex51_config.csv", "r"))


no_at_nums = ['18', '19', '20', '21', '22', '23', '32', '33', '34', '35', '36', '37']

# moving no_at files (empty multiple cloning site) to own directory
# for hiseq run 1
no_at_fs_m1 = [f for f in listdir(T_SINGLE_2_FILTERED_DIR_M1) if f[len('190716Lau_D19-79'): len('190716Lau_D19-79')+2] in no_at_nums] 

for f in no_at_fs_m1:
    mv_cmd = 'mv {}{} {}{}'.format(T_SINGLE_2_FILTERED_DIR_M1, f, T_SINGLE_2_FILTERED_DIR_NO_AT_M1,f)
    os.system(mv_cmd)
# for hiseq run 2 
no_at_fs_m2 = [f for f in listdir(T_SINGLE_2_FILTERED_DIR_M2) if f[len('190716Lau_D19-79'): len('190716Lau_D19-79')+2] in no_at_nums] 

for f in no_at_fs_m2:
    mv_cmd = 'mv {}{} {}{}'.format(T_SINGLE_2_FILTERED_DIR_M2, f, T_SINGLE_2_FILTERED_DIR_NO_AT_M2,f)
    os.system(mv_cmd)

# classifying files
class_samples(T_SINGLE_2_FILTERED_DIR_NO_AT_M1, T_SINGLE_2_CALLED_DIR_NO_AT_M1, df_config)
class_samples(T_SINGLE_2_FILTERED_DIR_NO_AT_M2, T_SINGLE_2_CALLED_DIR_NO_AT_M2, df_config)


