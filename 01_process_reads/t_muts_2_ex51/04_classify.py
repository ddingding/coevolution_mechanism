# classify reads and call single mutants

import yaml
import pandas as pd
from pipelineTools import class_samples
from constants import T_SINGLE_2_AT_SPLIT_DIR_M1, T_SINGLE_2_AT_SPLIT_DIR_M2, T_SINGLE_2_CALLED_DIR_M1, T_SINGLE_2_CALLED_DIR_M2

# read config files
config_dics = yaml.safe_load(open("./ex51_config.yaml", "r"))
df_config = pd.read_csv(open("./ex51_config.csv", "r"))


if __name__ == "__main__":
    # for hiseq run 1
    din = T_SINGLE_2_AT_SPLIT_DIR_M1
    dout = T_SINGLE_2_CALLED_DIR_M1
    class_samples(din, dout, df_config)

    # for hiseq run 2
    din = T_SINGLE_2_AT_SPLIT_DIR_M2
    dout = T_SINGLE_2_CALLED_DIR_M2
    class_samples(din, dout, df_config)
