# classify the reads into categories, such as single toxin variants, or multiple variants, read too short, etc.
import pandas as pd
from pipelineTools import class_samples
from constants import T_SINGLE_1_AT_SPLIT_DIR_M1, T_SINGLE_1_AT_SPLIT_DIR_M2, T_SINGLE_1_CALLED_DIR_M1, T_SINGLE_1_CALLED_DIR_M2
# read config files
df_config = pd.read_csv(open("./ex47_config.csv", "r"))

# for hiseq run 1
din = T_SINGLE_1_AT_SPLIT_DIR_M1
dout = T_SINGLE_1_CALLED_DIR_M1
class_samples(din, dout, df_config)

# for hiseq run 2
din = T_SINGLE_1_AT_SPLIT_DIR_M2
dout = T_SINGLE_1_CALLED_DIR_M2
class_samples(din, dout, df_config)
