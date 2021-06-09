# make sample objects and timecourse objects out of classified reads.
import libClass as lc
from os import listdir
from os.path import isfile, join

from pipelineTools import (
    map_primer_to_template,
    make_all_sample_obj,
    make_all_possible_tcs_single,
    pair_replicate_tc,
    pair_replicate_tc_single,
    load_all_pickles,
)

import plottingTools as pt
import libClassTools as lct

import yaml
import pandas as pd
import pickle


if __name__ == "__main__":
    # read config files
    config_dics = yaml.safe_load(open("./ex51_config.yaml", "r"))
    # py2
    # df_config = pd.read_csv(open('./ex51_config.csv', 'r'))

    # py3
    df_config = pd.read_csv(open("./ex51_config.csv", "r", encoding="ISO-8859-1"))

    # class_din = '/n/groups/marks/users/david/ex51/05_class_concat/'
    class_din = "/n/groups/marks/users/david/ex51/05_class_concat/mcs_renamed/"
    pickle_out = "/n/groups/marks/users/david/ex51/06_pickles/"
    plot_out = "/n/groups/marks/users/david/ex51/07_plots/"

    """
    #making samples objects. had to do this 4 times to get through all the sample opjects, memory issues.
    dic_samples = make_all_sample_obj(class_din, pickle_out, df_config)
    """

    """
    #and loading all the pickled sample objects in, doesn't work.
    print('trying to load all samples into mem...')
    sample_list = load_all_pickles(pickle_out, endswith_str = 'sample_obj.p')
    print('...done.')
    """

    s_name_list = [
        f.rstrip("_sample_obj.p")
        for f in listdir(pickle_out)
        if f.endswith("_sample_obj.p")
    ]

    # making all the tc objects and pickling
    print("s_name_list", sorted(s_name_list))
    tc_list = make_all_possible_tcs_single(s_name_list, pickle_out, df_config)
