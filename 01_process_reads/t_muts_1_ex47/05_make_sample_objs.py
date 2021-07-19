# take classified read output and create Sample objects, and then Timecourse objects

from pipelineTools import (
    map_primer_to_template,
    make_all_sample_obj,
    make_all_possible_tcs,
    pair_replicate_tc,
)

import yaml
import pandas as pd


# read config files
config_dics = yaml.safe_load(open("./ex47_config.yaml", "r"))
df_config = pd.read_csv(open("./ex47_config.csv", "r"))

class_din = "/n/groups/marks/users/david/ex47/05_class_concat/"
pickle_out = "/n/groups/marks/users/david/ex47/06_pickles/pickles_concat/"
plot_out = "/n/groups/marks/users/david/ex47/07_plots/concat/"

# making samples objects.
dic_samples = make_all_sample_obj(class_din, pickle_out, df_config)
# then make all the possible timecourses
tc_list = make_all_possible_tcs(dic_samples, pickle_out, df_config)
sample_list = list(dic_samples.values())

#plot some sequencing qc numbers
for s_obj in sample_list:
    pt.plot_sm_counts(s_obj, dir_out = plot_out)
    pt.plotSequencingQC(s_obj, dir_out = plot_out)
