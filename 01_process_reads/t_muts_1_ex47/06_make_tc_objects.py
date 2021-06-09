import pandas as pd
import matplotlib.pyplot as plt
import pickle
import libClassTools as lct

#%load_ext autoreload
import importlib
import pipelineTools as pit
import plottingTools as pt
from copy import deepcopy
from os.path import isfile
from os import listdir
from os.path import isfile, join


plot_out = "/n/groups/marks/users/david/ex47/07_plots/concat/"
pickleDir = "/n/groups/marks/users/david/ex47/06_pickles/pickles_concat/"
df_config = pd.read_csv(open("./ex47_config.csv", "r"))
dout = plot_out

"""
tc_objs = pit.load_all_pickles(pickleDir, 'tc.p')

#normalizing to each samples stop codons

for tc in tc_objs:
    tc_name = '_'.join(tc.samples[0].sample_n.split('_')[:-1])

    if not isfile(pickleDir+tc_name+'_tc_stop_fit_each.p'):

        tc.calculateFitness(first_timepoint=tc.samples[0].timepoint, fit_wrt_sample_stop=True)

        pickle.dump(tc, open(pickleDir+tc_name+'_tc_stop_fit_each.p', 'wb'))

        pt.make_tc_plots(tc, dout)
"""
#
tcs = [
    f for f in listdir(pickleDir) if isfile(join(pickleDir, f)) and f.endswith("tc.p")
]
print(tcs)


def pair_tc_fnames(tc_f_list, df_config):
    # expect list of filenames , like 179_W59T_tc_stop_fit_each.p
    # and returns a list of tuples of filenames for the replicates.
    df_rep = df_config[["primer", "replicate"]].dropna()

    tc_paired_list = []
    for tc_f in tc_f_list:
        pri1 = int(tc_f[:3])
        # check if the primer number is in primer (the columns are not symmetric)
        if pri1 in list(df_rep["primer"]):
            pri_rep = str(
                int(df_rep.loc[df_rep["primer"] == int(pri1), "replicate"].iloc[0])
            )
            tc_paired_list.append([tc_f, pri_rep + tc_f[3:]])
    return tc_paired_list


paired_tc_fnames = pair_tc_fnames(tcs, df_config)
print(paired_tc_fnames)

for [tc1_f, tc2_f] in paired_tc_fnames:
    print(tc1_f, tc2_f)

    tc1 = pickle.load(open(pickleDir + tc1_f, "rb"))
    tc2 = pickle.load(open(pickleDir + tc2_f, "rb"))

    # make the heatmaps
    pt.make_tc_plots(tc1, plot_out)
    pt.make_tc_plots(tc2, plot_out)

    # get the fitness correlation between replicates
    tc1_t = int(df_config.loc[df_config.primer == int(tc1.samples[-1].sample_n[:3])].t)
    tc2_t = int(df_config.loc[df_config.primer == int(tc2.samples[-1].sample_n[:3])].t)
    pt.get_fit_corr_tcs(
        tc1,
        tc2,
        t=tc1_t,
        outdir=plot_out + "rep/",
        plot_aa_fitness=True,
        figname="aa_corr " + tc1.samples[-1].sample_n,
    )
    pt.get_fit_corr_tcs(
        tc1,
        tc2,
        t=tc1_t,
        outdir=plot_out + "rep/",
        plot_aa_fitness=False,
        figname="corr " + tc1.samples[-1].sample_n,
    )

    ave_sample = lct.get_ave_fit_sample([tc1.samples[-1], tc2.samples[-1]])
    pickle.dump(
        ave_sample, open(pickleDir + ave_sample.sample_n + "_sample_ave.p", "wb")
    )

    # plot average fitness plots
    pt.plot_fitness_hist(ave_sample, dir_out=dout + "ave/")
    pt.plot_fitness_heatmap(ave_sample, dir_out=dout + "ave/")

"""
tc_objs = pit.load_all_pickles(pickleDir, 'tc_stop_fit_each.p')

# making the average fitness samples.
rep_tcs = pit.pair_replicate_tc(tc_objs, df_config)
ave_sample_dic = {}
for (tc1, tc2) in rep_tcs:
    ave_sample = lct.get_ave_fit_sample([tc1.samples[-1], tc2.samples[-1]])
    pickle.dump(ave_sample, open(pickleDir+ave_sample.sample_n+'.p', 'wb'))
    ave_sample_dic[ave_sample.sample_n] = ave_sample
"""
