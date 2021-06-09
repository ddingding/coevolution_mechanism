import libClass as lc
from os import listdir
from os.path import isfile, join

from pipelineTools import (
    map_primer_to_template,
    make_all_sample_obj,
    make_all_possible_tcs,
    pair_replicate_tc,
)

import plottingTools as pt

import yaml
import pandas as pd
import pickle

# Goal:
# 1 make sample objects.
# save as pickles


# plots


if __name__ == "__main__":
    # read config files
    config_dics = yaml.safe_load(open("./ex47_config.yaml", "r"))
    df_config = pd.read_csv(open("./ex47_config.csv", "r"))
    class_din = "/n/groups/marks/users/david/ex47/05_class_concat/"
    pickle_out = "/n/groups/marks/users/david/ex47/06_pickles/pickles_concat/"
    plot_out = "/n/groups/marks/users/david/ex47/07_plots/concat/"

    # making samples objects.
    """
    dic_samples = make_all_sample_obj(class_din, pickle_out, df_config)
    tc_list = make_all_possible_tcs(dic_samples, pickle_out, df_config)
    sample_list = list(dic_samples.values())
    #sample_list = load_all_pickles(pickle_out, '178_E79H_sample_obj.p')

    for s_obj in sample_list:
        pt.plot_sm_counts(s_obj, dir_out = plot_out)
        pt.plotSequencingQC(s_obj, dir_out = plot_out)
    """
    # making timecourses
    tc_list = load_all_pickles(pickle_out, "_tc.p")

    # make all the tc plots

    for tc in tc_list:
        pt.make_tc_plots(tc, plot_out)

    # import pdb; pdb.set_trace()
    ##
    paired_tc_list = pair_replicate_tc(tc_list, df_config)
    for (tc1, tc2) in paired_tc_list:

        tc1_t = int(
            df_config.loc[df_config.primer == int(tc1.samples[-1].sample_n[:3])].t
        )
        tc2_t = int(
            df_config.loc[df_config.primer == int(tc2.samples[-1].sample_n[:3])].t
        )
        assert tc1_t == tc2_t

        pt.getFitCorr(
            tc1,
            tc2,
            t=tc1_t,
            outdir=plot_out + "rep/",
            plotAAFitness=True,
            figname="aa_corr " + tc1.samples[-1].sample_n,
        )
        pt.getFitCorr(
            tc1,
            tc2,
            t=tc1_t,
            outdir=plot_out + "rep/",
            plotAAFitness=False,
            figname="corr " + tc1.samples[-1].sample_n,
        )
