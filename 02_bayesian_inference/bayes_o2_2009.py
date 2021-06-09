'''
This script takes MCMC arguments including filenames for the raw data, and performs MCMC inference using PyStan to
infer growth rates for amino acid mutants.
'''

import sys
import stanTools
import stanModels as smods
import pystan
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
from os.path import isfile
import numpy as np
import csv
import pandas as pd

if len(sys.argv) != 9:
    print('too many or few args given to bayes_o2.py: {}'.format(' '.join(map(str, sys.argv))))

# reading MCMC parameters from args
pickle_dir = sys.argv[3]
plot_out = sys.argv[4]
mcmc_iter_all = int(sys.argv[5])
rand_seed = int(sys.argv[6])
n_chains = int(sys.argv[7])
init_r = float(sys.argv[8])
mcmc_n = n_chains
thin = 2

# fetching experiment raw count data
f1 = sys.argv[1]
f2 = sys.argv[2]
data_all, data_r1, data_r2, aa_dic, df_m, df_syn = stanTools.get_stan_data(f1, f2)

f_name1 = f1.split('/')[-1][:-4]  # get rid of '.csv'
f_name2 = f2.split('/')[-1][:-4]
f_name_rep = f_name1[:-3]  # get rid of '_tc' for the sm8 files

# aa_dic is a dictionary of integer to a particular amino acid mutant. This ordering is used in stan
# write out the aa dic for later reference
fout = pickle_dir + "aa_dic" + f_name_rep + ".csv"
with open(fout, "w") as fout:
    w = csv.writer(fout)
    for key, val in aa_dic.items():
        w.writerow([key, val])

# load and compile the the stan model
sm_file = pickle_dir + 'sm_bhs_' + f_name_rep + '.pkl'
if not isfile(sm_file):
    sm = pystan.StanModel(model_code=smods.bayes_hierarchical_syn)
    with open(sm_file, 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(sm_file, 'rb'))

# sampling
# had to set n_jobs=1, meaning sequential chain sampling, due to some Multi-thread/pooling issues.
f_file = pickle_dir + 'f_sm_bhs_' + f_name_rep + '.p'
if not isfile(f_file):
    if init_r == None:
        f_sm = sm.sampling(data=data_all, iter=mcmc_iter_all, chains=mcmc_n, seed=rand_seed, n_jobs=1,
                           thin=thin, control={'max_treedepth': 15, 'adapt_delta': 0.99})
    else:
        f_sm = sm.sampling(data=data_all, iter=mcmc_iter_all, chains=mcmc_n, seed=rand_seed, init_r=init_r,
                           n_jobs=1,
                           thin=thin, control={'max_treedepth': 15, 'adapt_delta': 0.99})
    with open(f_file, 'wb') as f:
        pickle.dump(f_sm, f)
else:
    f_sm = pickle.load(open(f_file, 'rb'))

# create an output file. Note that all the results in this file are in terms of ln with base e, so need to convert to
# log_2
df_fit_rep, df_fit_mean_rep = stanTools.dump_and_convert(f_sm, 'sm_bhs_diff_w_aa_m_' + f_name_rep, aa_dic,
                                                         pickle_dir,
                                                         col_pre='diff_w_aa_m')

##### to plot some posterior predictive checks
# get real data
print('fetching observed data....')
data_all, aa_dic, df_m, df_syn = stanTools.get_formatted_stan_data(f1, f2)  # df_syn needs a read ratio

# get the inferred df
print('loading dataframe calculated....')
f_name_load = 'df_sm_bhs_diff_w_aa_m_' + f_name_rep + '.csv'
df_fit_rep = pd.read_csv(pickle_dir + f_name_load)

# add lrr_syn1', 'lrr_syn2', 'lrr_m1', 'lrr_m2', 'lrr_aa_m1', 'lrr_aa_m2', 'lrr_aa_m_12' to data_subset, data_all
print('adding lrrs to observed data...')
data_all = stanTools.add_codon_lrr_to_stan_data(data_all)
data_all = stanTools.add_observed_aa_lrr_to_stan_data(df_m, data_all)

print('adding lrrs to inferred data...')
df_fit_rep = stanTools.add_posterior_codon_lrr_to_df_fit(df_fit_rep, data_all)
df_fit_rep = stanTools.add_posterior_aa_lrr(df_fit_rep, data_all)
# plotting PPCs
fout = pickle_dir + 'ppc_sm_bhs_' + f_name_rep + '.svg'
stanTools.ppc_figures_one_df_fit(df_fit_rep, data_all, df_m, df_syn, fout=fout)
