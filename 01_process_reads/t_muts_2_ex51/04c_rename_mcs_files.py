import libClass as lc
from os import listdir
from os.path import isfile, join
import os

from pipelineTools import (
    map_primer_to_template,
    make_all_sample_obj,
    make_all_possible_tcs_single,
    pair_replicate_tc,
    pair_replicate_tc_single,
    load_all_pickles,
)

import plottingTools as pt

import yaml
import pandas as pd
import pickle
import numpy as np


mcs_class_in = "/n/groups/marks/users/david/ex51/05_class_concat/mcs/"
mcs_class_rename_out = "/n/groups/marks/users/david/ex51/05_class_concat/mcs_renamed/"

fs = [f for f in listdir(mcs_class_in)]


df_config = pd.read_csv(open("./ex51_config.csv", "r", encoding="ISO-8859-1"))
bmc_to_pri = dict(zip(df_config["bmc_index_lane_1"], df_config["primer"]))

print(bmc_to_pri)


lane_to_pri = {}
for bmc, pri in bmc_to_pri.items():
    if not np.isnan(pri):
        l = bmc.split(" : ")[0][-10:-2]
        lane_to_pri[l] = str(int(pri))
print("lane_to_pri", lane_to_pri, fs)

for f in fs:
    print(f)
    l = f.split("_")[1][:8]
    pri = lane_to_pri[l]
    cmd = "cp {}{} {}{}_mcs_class.tsv".format(
        mcs_class_in, f, mcs_class_rename_out, pri
    )
    print(cmd)
    os.system(cmd)
