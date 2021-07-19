'''
tools for dealing with stan and the models

functions for
1) formatting read data into dictionaries that can be fed into the stan model for inference
2) analysing/summarizing stan inferred data
3) calculating posterior predictive checks
4) simulating synthetic datasets to check that the stan model infers values that make sense

'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
import scipy.stats as st

import libClassTools as lct


#########################################################################
# formatting data into a dictionary that is compatible with feeding into stan model 'bayes_hierarchical_syn'
def get_stan_data(pin, pin2, pseudo_count=0):
    '''

    :param pin: file path for the first replicate
    :param pin2: filepath for the second replicate
    :param pseudo_count:
    :return:
    '''
    # merge dfs
    df_reps = lct.merge_df_counts(pin, pin2)

    # dropping nan in raw_count_r1, and fill the nans in raw_count_post_r1 (these are non detected in post selection)
    df_rep1 = df_reps[['mutkey', 'codonkey', 'raw_count_r1', 'raw_count_post_r1']]
    df_rep2 = df_reps[['mutkey', 'codonkey', 'raw_count_r2', 'raw_count_post_r2']]

    df_rep1['raw_count_post_r1'] = df_rep1['raw_count_post_r1'].fillna(0)
    df_rep2['raw_count_post_r2'] = df_rep2['raw_count_post_r2'].fillna(0)
    df_rep1 = df_rep1.dropna()
    df_rep2 = df_rep2.dropna()
    df_reps = df_rep1.merge(df_rep2, left_on='codonkey', right_on='codonkey', suffixes=('', '_x'))

    # sorting by mutkey
    df_count = df_reps[['mutkey', 'codonkey', 'raw_count_r1', 'raw_count_post_r1', 'raw_count_r2', 'raw_count_post_r2']]
    df_count = df_count.sort_values(by='mutkey')

    # split into synonymous wt mutants and mutant counts
    df_m = df_count.loc[df_count["mutkey"].str[0] != df_count["mutkey"].str[-1]]
    df_syn = df_count.loc[df_count["mutkey"].str[0] == df_count["mutkey"].str[-1]]
    d, df_m = get_mut_counts(df_m)

    # adding pseudocount
    df_m[['raw_count_r1', 'raw_count_post_r1', 'raw_count_r2', 'raw_count_post_r2']] = df_m[['raw_count_r1',
                                                                                             'raw_count_post_r1',
                                                                                             'raw_count_r2',
                                                                                             'raw_count_post_r2']] + \
                                                                                       pseudo_count
    df_syn[['raw_count_r1', 'raw_count_post_r1', 'raw_count_r2', 'raw_count_post_r2']] = df_syn[['raw_count_r1',
                                                                                                 'raw_count_post_r1',
                                                                                                 'raw_count_r2',
                                                                                                 'raw_count_post_r2']] + \
                                                                                         pseudo_count

    data_all = {
        'K'       : int(len(df_syn.raw_count_r1)),
        'c_pre1'  : df_syn.raw_count_r1.astype(int),
        'c_aft1'  : df_syn.raw_count_post_r1.astype(int),
        'c_pre2'  : df_syn.raw_count_r2.astype(int),
        'c_aft2'  : df_syn.raw_count_post_r2.astype(int),
        'Km'      : int(len(d)),
        'N'       : len(df_m),
        's'       : df_m.g_idx,
        'c_pre_m1': df_m.raw_count_r1.astype(int),
        'c_aft_m1': df_m.raw_count_post_r1.astype(int),
        'c_pre_m2': df_m.raw_count_r2.astype(int),
        'c_aft_m2': df_m.raw_count_post_r2.astype(int)

    }

    data_rep_1 = {'K'      : int(len(df_syn.raw_count_r1)),
                  'c_pre'  : df_syn.raw_count_r1.astype(int),
                  'c_aft'  : df_syn.raw_count_post_r1.astype(int),
                  'Km'     : int(len(d)),
                  'N'      : len(df_m),
                  'c_pre_m': df_m.raw_count_r1.astype(int),
                  'c_aft_m': df_m.raw_count_post_r1.astype(int),
                  's'      : df_m.g_idx
                  }

    data_rep_2 = {
        'K'      : int(len(df_syn.raw_count_r2)),
        'c_pre'  : df_syn.raw_count_r2.astype(int),
        'c_aft'  : df_syn.raw_count_post_r2.astype(int),
        'Km'     : int(len(d)),
        'N'      : len(df_m),
        'c_pre_m': df_m.raw_count_r2.astype(int),
        'c_aft_m': df_m.raw_count_post_r2.astype(int),
        's'      : df_m.g_idx

    }
    return data_all, data_rep_1, data_rep_2, d, df_m, df_syn


def get_formatted_stan_data(f1, f2):
    # calls get_stan_data, but also adds adds log read ratios,
    # scales lrr2 by the sequencing depth of the different replicates

    # 1: fetch data
    data_all, data_r1, data_r2, aa_dic, df_m, df_syn = get_stan_data(f1, f2, pseudo_count=0)

    # 2: add relevant columns
    ## creating df_m with all the read ratios and variances
    # rename columns so that they're compatible with add_ratios_lrr_df_m
    df_m = df_m.rename(columns={
        'raw_count_r1'     : 'c_pre_m1',
        'raw_count_post_r1': 'c_aft_m1',
        'raw_count_r2'     : 'c_pre_m2',
        'raw_count_post_r2': 'c_aft_m2'
    })

    df_m = add_ratios_lrr_df_m(df_m)

    # to get lrr2 scaled
    # to get the ratios:
    df_sum = np.sum(df_m)
    read_ratios = (df_sum['c_aft_m1'] / df_sum['c_pre_m1']) / (df_sum['c_aft_m2'] / df_sum['c_pre_m2'])

    df_m['lrr2_scaled'] = df_m['lrr2'] - np.log(1 / read_ratios)

    '''
    # to verify the scaling worked in the right direction
    plt.figure()
    plt.scatter(df_m.lrr1, df_m.lrr2_scaled)
    plt.plot([-5,1],[-5,1], c='black')
    plt.show()
    '''

    # 3 adding synonymous growth rates to the

    df_syn = df_syn.rename(columns={'raw_count_r1'     : 'c_pre1',
                                    'raw_count_r2'     : 'c_pre2',
                                    'raw_count_post_r1': 'c_aft1',
                                    'raw_count_post_r2': 'c_aft2'})

    df_syn['ratio1'] = df_syn['c_aft1'] / df_syn['c_pre1']
    df_syn['ratio2'] = df_syn['c_aft2'] / df_syn['c_pre2']  # found the bug!
    df_syn['lrr1'] = np.log(df_syn['ratio1'])
    df_syn['lrr2'] = np.log(df_syn['ratio2'])

    # just to check the overall shape of the distribution, the values don't mean much otherwise.
    # to viz the distribution from which I'm sampling growth ratios
    plt.figure()
    plt.hist(np.ma.masked_invalid(df_m.lrr1))
    plt.hist(np.ma.masked_invalid(df_m.lrr2))

    plt.hist(np.ma.masked_invalid(df_syn.lrr1))

    plt.show()

    return data_all, aa_dic, df_m, df_syn


# todo check that the g_idx are actually in ascending incremental order in the resulting dataframe. ie. the sort here
#  is not different than the df_m sort by mutkey
def get_mut_counts(df_m):
    # stan needs to know which codon mutations correspond to a particular amino acid substitution to build a
    # hierarchy on
    # this dictionary d maps the the amino acid mutation ('mutkey') to an integer, and adds this integer to the
    # mutation dataframe df_m
    d = dict(zip(sorted(set(df_m['mutkey'])), np.array(range(len(set(df_m['mutkey'])))) + 1))

    # create index column for stan, each integer coding for a particular mutant.
    df_m['g_idx'] = df_m['mutkey'].map(d)

    return d, df_m


def dump_and_convert(fit, save_name, aa_dic, pickle_dir, col_pre='diff_r', ):
    '''
    take a stan fit object
        - save it as pickle
        - convert to dataframe, and save that
        - calculate the mean for columns starting with 'col_pre', which are the codon mutant or amino acid mutant
        columns
    :param fit: stan fit object
    :param save_name: should be whatever comes after f_, and df_ for saving
    :param aa_dic: the dictionary from get_mut_counts() to allow referencing amino acid mutant index and the mutant
    :param col_pre: stan output columns prefic. stan spits out columns named diff_r[n], where n is the amino acid
            numbering from aa_dic.
    :param pickle_dir: directory to save the dataframe to

    '''
    print('dumping ' + save_name + '...')

    # saving and converting to pickles
    # pickle.dump(fit, open(pickle_dir + 'f_' + save_name + '.p', 'wb'))
    df_fit = fit.to_dataframe()
    print('saving.')
    df_fit.to_csv(pickle_dir + 'df_' + save_name + '.csv')

    # saving the dataframe of interest
    df_mean_fit = create_df_mean_hdi(df_fit, aa_dic, col_pre=col_pre)
    df_mean_fit.to_csv(pickle_dir + 'df_mean' + save_name + '.csv')
    return df_fit, df_mean_fit


def create_df_mean_hdi(df_fit, dic_aakey, col_pre='diff_r', exp=False):
    '''
    takes a df_fit from stan (containing the psoterior values of all MCMC chain samples),
     and returns a dataframe of the mean and 95% Highest density interval of the posterior for each growth rate
    :param dic_aakey: the dictionary from get_mut_counts() to allow referencing amino acid mutant index and the mutant
    returns df with columns mean_fit, hdi_lower, hdi_upper
    '''
    df = create_aamut_col_df(dic_aakey, df_fit, col_pre=col_pre)
    samples = df_fit[df.col.values]
    if exp:
        samples = np.exp(samples)
    df['mean_fit'] = np.mean(samples, axis=0).values

    df['hdi_lower'] = (samples.apply(compute_HDI_lower, axis=0, raw=True)).values
    df['hdi_upper'] = (samples.apply(compute_HDI_upper, axis=0, raw=True)).values
    return df


def create_aamut_col_df(d, df_fit, col_pre='diff_r'):
    ''' create df with aamut and column name to index aamut
    :param d:       the dictionary from get_mut_counts() to allow referencing amino acid mutant index and the mutant
    :param df_fit:  the df_fit dataframe from stan
    '''
    df_mut = pd.DataFrame.from_dict(list(d.items()), orient='columns')
    df_mut = df_mut.rename(columns={0: 'aa_mut', 1: 'num'})

    # get diff_r columns
    diff_cols = [c for c in df_fit.columns if c.startswith(
        col_pre + '[')]  # the bracket search so that I can search for w_m[1], and exclude w_m1 and w_m1_std
    # make dictionary of aamut number to diff_r
    num_to_col = dict(zip(map(int, [c.lstrip(col_pre)[1:-1] for c in diff_cols]), diff_cols))
    # convert to dataframe
    df_col = pd.DataFrame.from_dict(list(num_to_col.items()))
    df_col = df_col.rename(columns={0: 'num', 1: 'col'})
    df_merge = df_mut.merge(df_col, left_on='num', right_on='num')
    return df_merge


def compute_HDI_upper(chain, interval=.95):
    '''wrapper for compute_HDI'''
    return compute_HDI(chain, interval=interval)['Upper']


def compute_HDI_lower(chain, interval=.95):
    '''wrapper for compute_HDI'''
    return compute_HDI(chain, interval=interval)['Lower']


def compute_HDI(chain, interval=.95):
    '''
    computes highest density interval, HDI,  that spans 95% of MCMC samples for a particular amino acid or codon mutant
    '''
    # sort chain using the first axis which is the chain
    chain.sort()
    # how many samples did you generate?
    nSample = chain.size
    # how many samples must go in the HDI?
    nSampleCred = int(np.ceil(nSample * interval))
    # number of intervals to be compared
    nCI = nSample - nSampleCred
    # width of every proposed interval
    width = np.array([chain[i + nSampleCred] - chain[i] for i in range(nCI)])
    # index of lower bound of shortest interval (which is the HDI)
    best = width.argmin()
    # put it in a dictionary
    HDI = {'Lower': chain[best], 'Upper': chain[best + nSampleCred], 'Width': width.min()}
    return HDI


def calc_p_val_from_threshold(chain, threshold=0):
    ''' takes a chain and calculates the p-vale of being at the threshold 0
    default threshold is 0, which is where the growth rate of an antitoxin single mutant with a wild-type toxin falls
    '''
    # print(chain[:5])
    chain_no_nan = chain.dropna()
    # print(len(chain), len(chain_no_nan))
    total_n = len(chain_no_nan)

    # check there are enough samples to still calculate a good score
    if total_n < len(chain) / 10:
        return np.nan

    # check they are not all just +inf values, ie. occurs when observed post selection read is 0.
    if np.sum(chain_no_nan == np.inf) == total_n or np.sum(chain_no_nan == -np.inf) == total_n:
        return np.nan

    vals_below = np.sum(chain_no_nan.values < threshold)

    if vals_below == 0:
        return 0
        # return (1 / total_n)

    elif vals_below == total_n:
        return 1
        # return ((vals_below - 1) / total_n)

    else:
        return vals_below / total_n


def create_df_mean_percentiles(df_fit, dic_aakey, col_pre='diff_w_aa_m'):
    '''
    function to calculate [mean, 2.5%, 5%, 10, 16, 25, 30, 50, 70, 75, 84, 90, 95, 97.5] percentiles from df_fit
    '''
    df = create_aamut_col_df(dic_aakey, df_fit, col_pre=col_pre)
    samples = df_fit[df.col.values]

    df['mean_fit'] = np.mean(samples, axis=0).values

    percents = [2.5, 5, 10, 16, 25, 30, 50, 70, 75, 84, 90, 95, 97.5]
    for perc in percents:
        df[str(perc) + '%'] = np.percentile(samples, perc, axis=0)
    p_vals = samples.apply(calc_p_val_from_threshold, axis=0).values
    df['p'] = p_vals
    return df


def create_df_perc_epi_p_log2(f, pin, fout_name, sample_str, df_epi, col_pre='diff_w_aa_m'):
    '''
    load stan df_fit (containing all the mcmc samples, in terms of ln()) and convert to df percentage, which is a
    summary of all the
    percentiles, where delta_growth_rates are measured in log2()
    :param f:
    :param pin:
    :param fout_name:
    :param sample_str:
    :param df_epi: contains a column 'yhat_double_minus_yobs_single' that indicates the indpendent, nonlinear double
    mutant expectation of the growth rate
    :return: df_percentage
    # includes a p value for being epistatic or not, depending on how many of the samples fall above or below the
    # independent prediction
    '''
    # load aa dic
    print('loading aa_dic')
    aa_f = 'aa_dic' + sample_str + '.csv'
    df_aa_dic = pd.read_csv(pin + aa_f, header=None)
    aa_dic = dict(zip(df_aa_dic[0], df_aa_dic[1]))

    # load df
    print('loading df')
    df_fit = pd.read_csv(pin + f)

    # get a df of t_mut to col_name: like
    '''aa_mut	num	col
        0	A10C	1	diff_w_aa_m[1]
        1	A10D	2	diff_w_aa_m[2]
    '''
    df = create_aamut_col_df(aa_dic, df_fit, col_pre=col_pre)

    diff_w_cols = df.col.values
    # convert to log 2 base, not natural base
    df_fit[diff_w_cols] = df_fit[diff_w_cols] / np.log(2)

    # get the values
    samples = df_fit[diff_w_cols]

    df['mean_fit'] = np.mean(samples, axis=0).values
    percents = [2.5, 5, 10, 16, 25, 30, 50, 70, 75, 84, 90, 95, 97.5]
    for perc in percents:
        df[str(perc) + '%'] = np.percentile(samples, perc, axis=0)

    # calculate p-value from low values.
    print('calculating p value from 0 (AT single)')
    # calculate p_value from 0 (at single)
    p_vals = samples.apply(calc_p_val_from_threshold, axis=0).values
    df['p'] = p_vals

    ########## calculate p_vals for independent
    print('calculating epistatic p-value')
    # get the antitoxin specific double mutants
    at_mut = sample_str.split('_')[-1]

    df_epi_at = df_epi.loc[df_epi.at_mut == at_mut]

    yobs_single = df_epi.loc[(df_epi.t_mut == 'wtT') & (df_epi.at_mut == at_mut)].yobs_from_wt.values[0]
    print('yobs_single', yobs_single)
    # residual of independent double mutant prediction(yhat_obs_from_wt ) minus what the observed yobs_single AT
    # mutant is
    df_epi_at['yhat_double_minus_yobs_single'] = np.subtract(df_epi_at.yhat_from_wt.values, yobs_single)

    # get the residuals in the correct order as samples
    yhat_minus_yobs = []
    for aa_mut in df.aa_mut:
        # check if residual values are missing or not. if so, just add a nan value.
        if aa_mut in df_epi_at.t_mut.values:
            res = df_epi_at.loc[df_epi_at.t_mut == aa_mut].yhat_double_minus_yobs_single.values[0]
            yhat_minus_yobs.append(res)
        else:
            yhat_minus_yobs.append(np.nan)
    yhat_minus_yobs = np.array(yhat_minus_yobs)

    # subtracting from samples what each columns difference, then I can just count in each column which mutants are
    # above 0, or below to get a p-value
    samples_from_indep = np.subtract(samples, yhat_minus_yobs)
    p_vals_epi = samples_from_indep.apply(calc_p_val_from_threshold, axis=0).values
    df['p_epi'] = p_vals_epi

    # print(p_vals)
    print('writing rescaled_mean')
    df.to_csv(pin + fout_name)
    return df, df_fit


######################################################################
# posterior predictive checking functions
##########################################################################
def add_codon_lrr_to_stan_data(stan_data):
    '''adds log read ratios to the observed stan data dictionary'''
    # true log read ratios of all the synonymous mutants
    stan_data['lrr_syn1'] = np.log(stan_data['c_aft1'] / stan_data['c_pre1'])
    stan_data['lrr_syn2'] = np.log(stan_data['c_aft2'] / stan_data['c_pre2'])

    # log read ratios of all the codons
    stan_data['lrr_m1'] = np.log(stan_data['c_aft_m1'] / stan_data['c_pre_m1'])
    stan_data['lrr_m2'] = np.log(stan_data['c_aft_m2'] / stan_data['c_pre_m2'])

    return stan_data


def add_observed_aa_lrr_to_stan_data(df_m, data_dic):
    '''
    adding the observed df_m in the dataframe version of stan_data dictionary (data_dic)

    '''
    df_aa_mean = df_m.groupby('g_idx').mean()
    data_dic['lrr_aa_m1'] = df_aa_mean.lrr1
    data_dic['lrr_aa_m2'] = df_aa_mean.lrr2
    data_dic['lrr_aa_m_12'] = np.mean(df_aa_mean[['lrr1', 'lrr2']], axis=1)
    return data_dic


def add_posterior_codon_lrr_to_df_fit(df_fit, stan_data):
    # adding log read ratio columns to the generated df_m
    # generated quantities are c_pre1_pred, c_aft1_pred, c_pre_m1_pred, c_aft_m1_pred

    df_m_curr = df_fit

    input_data = stan_data
    n_syn = input_data['K']
    n_codon = input_data['N']

    # for synonymous mutants
    for i in range(1, n_syn + 1):
        df_m_curr['lrr_syn1_pred[{}]'.format(i)] = np.log(
            df_m_curr['c_aft1_pred[{}]'.format(i)] / df_m_curr['c_pre1_pred[{}]'.format(i)])
        df_m_curr['lrr_syn2_pred[{}]'.format(i)] = np.log(
            df_m_curr['c_aft2_pred[{}]'.format(i)] / df_m_curr['c_pre2_pred[{}]'.format(i)])

    for i in range(1, n_codon + 1):
        df_m_curr['lrr_m1_pred[{}]'.format(i)] = np.log(
            df_m_curr['c_aft_m1_pred[{}]'.format(i)] / df_m_curr['c_pre_m1_pred[{}]'.format(i)])
        df_m_curr['lrr_m2_pred[{}]'.format(i)] = np.log(
            df_m_curr['c_aft_m2_pred[{}]'.format(i)] / df_m_curr['c_pre_m2_pred[{}]'.format(i)])

    return df_m_curr


def add_posterior_aa_lrr(df_fit, data):
    ''' adds amino acid average log read ratios per replicate, as well as across replicates.
    '''

    s = sorted(data['s'].values)
    n_codons = data['N']
    n_aa = data['Km']

    def get_aa_mean(df_fit, pref_cols_in, pref_cols_out):
        # select columns of interest in order
        cois = [pref_cols_in + '[{}]'.format(i) for i in range(1, n_codons + 1)]
        df_lrr_codon = df_fit[cois]
        # transpose to groupby s, and take mean
        df_lrr_codon_t = df_lrr_codon.T
        df_lrr_codon_t['s'] = s
        df_lrr_aa_mean = df_lrr_codon_t.groupby('s').mean()
        # change the index to the preferred coloumn names, and transpose
        df_lrr_aa_mean.index = [pref_cols_out + '[{}]'.format(i) for i in range(1, n_aa + 1)]
        df_lrr_aa_mean_t = df_lrr_aa_mean.T
        df_fit = pd.concat([df_fit, df_lrr_aa_mean_t], axis=1)
        return df_fit

    df_fit = get_aa_mean(df_fit, 'lrr_m1_pred', 'lrr_aa_m1_pred')
    df_fit = get_aa_mean(df_fit, 'lrr_m2_pred', 'lrr_aa_m2_pred')

    # get mean of amino acids
    mean_aa_dic = {}
    for i in range(1, n_aa + 1):
        mean_aa_dic['lrr_aa_m_12_pred[' + str(i) + ']'] = np.mean(
            df_fit[['lrr_aa_m1_pred[{}]'.format(i), 'lrr_aa_m2_pred[{}]'.format(i)]], axis=1)
    df_means = pd.DataFrame(mean_aa_dic)
    df_fit = pd.concat([df_fit, df_means], axis=1)

    return df_fit


def add_ratios_lrr_df_m(df_m,
                        col_pairs={'ratio1': ('c_aft_m1', 'c_pre_m1'),
                                   'ratio2': ('c_aft_m2', 'c_pre_m2')},
                        pseudocount=0,
                        correct_seq_depth=False
                        ):
    ''' take a dataframe df_m, and creates ratio and log ratio columns'''

    for new_col_n, (col_aft, col_pre) in col_pairs.items():
        if correct_seq_depth == False:
            df_m[new_col_n] = (df_m[col_aft] + pseudocount) / (df_m[col_pre] + pseudocount)
        else:
            df_m[new_col_n] = (df_m[col_aft] + pseudocount) / (df_m[col_pre] + pseudocount) * np.sum(
                df_m[col_pre] + pseudocount) / np.sum(df_m[col_aft] + pseudocount)

    df_m['lrr1'] = np.log(df_m['ratio1'])
    df_m['lrr2'] = np.log(df_m['ratio2'])

    return df_m


def ppc_figures_one_df_fit(df_fit, stan_data, df_stan_data, df_syn_stan_data, fout=None):
    '''

    plotting posterior predictive checks, ie. compare the distribution of simulated count data from the inferred
    parameters vs. the observed data to check that the model does not think the observed data is really unlikely
    plotting PPCs of each observed replicate dataset for
        wild-type synonymous log read ratios
        standard deviation of wild-type synonymous log read ratios
        codon mutant log read ratios
        standard deviations of synonymous amino acid mutant log read ratios
        mean synonymous amino acid mutant log read ratios

    :param df_fit: stan inferred data
    :param stan_data: dictionary of stan data
    :param df_stan_data: dataframe of stan data with log read ratios added
    :param df_syn_stan_data: dataframe of wild-type synonymous stan data with log read ratios added
    :param fout: output file to save plot
    '''

    n_rows = 7
    n_cols = 2
    fig = plt.figure(figsize=(n_cols * 3, n_rows * 2))
    gs = gridspec.GridSpec(n_rows, n_cols)
    plt.subplots_adjust(hspace=1, wspace=0.5)

    # rename for consistency
    df_syn_stan_data = df_syn_stan_data.rename(columns={'lrr1': 'lrr_syn1', 'lrr2': 'lrr_syn2'})
    df_stan_data = df_stan_data.rename(columns={'lrr1': 'lrr_m1', 'lrr2': 'lrr_m2'})

    # 1) synonymous lrrs biased?
    obs_vals = df_syn_stan_data.lrr_syn1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_syn1', ax=plt.subplot(gs[0]),
                               p_title='synonymous wt log read ratio ranks', title_size=8)

    obs_vals = df_syn_stan_data.lrr_syn2.values
    print('first five obs vals: lrr2_syn', obs_vals[:5])
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_syn2', ax=plt.subplot(gs[1]),
                               p_title='synonymous wt log read ratio ranks lrr2', title_size=8)

    # 2) is the standard dev of lrr of the synonymous fitting?
    plot_ppc(plt.subplot(gs[2]), df_fit, 'lrr_syn1', df_syn_stan_data, np.std,
             p_title='synonymous wt log read ratio stdev',
             xlabel='stdev')

    plot_ppc(plt.subplot(gs[3]), df_fit, 'lrr_syn2', df_syn_stan_data, np.std,
             p_title='synonymous wt log read ratio stdev',
             xlabel='stdev')

    # 3) mutant log read ratios have bias?
    obs_vals = df_stan_data.lrr_m1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_m1', ax=plt.subplot(gs[4]),
                               p_title='mutant log read ratio ranks', title_size=8)

    obs_vals = df_stan_data.lrr_m2.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_m2', ax=plt.subplot(gs[5]),
                               p_title='mutant log read ratio ranks', title_size=8)

    # 4) synonymous mutant standard deviations compared to observed ones
    obs_mutant_stds = df_stan_data.groupby('g_idx').std()
    df_cod_std_lrr1_diff_full = compare_posterior_std_synon(df_fit, stan_data, 'lrr_m1_pred[',
                                                            obs_mutant_stds['lrr_m1'], ax=plt.subplot(gs[6]),
                                                            p_title='synonymous mutant log read ratio stdev',
                                                            title_size=8)
    df_cod_std_lrr2_diff_full = compare_posterior_std_synon(df_fit, stan_data, 'lrr_m2_pred[',
                                                            obs_mutant_stds['lrr_m2'],
                                                            p_title='synonymous mutant log read ratio stdev',
                                                            title_size=8,
                                                            ax=plt.subplot(gs[7]))
    # 5) mean log read ratio of amino acid mutants
    obs_vals = stan_data['lrr_aa_m1'].values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_aa_m1', ax=plt.subplot(gs[8]),
                               p_title='mean amino acid mutant log read ratio')

    obs_vals = stan_data['lrr_aa_m2'].values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_aa_m2', ax=plt.subplot(gs[9]),
                               p_title='mean amino acid mutant log read ratio')

    ''' 
    #These are for counts
    #3) observed amino acid mutant counts
    obs_vals = df_stan_data.c_pre_m1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='c_pre_m1', ax=plt.subplot(gs[8]))

    obs_vals = df_stan_data.c_aft_m1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='c_aft_m1', ax=plt.subplot(gs[9]))

    obs_vals = df_stan_data.c_pre_m2.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='c_pre_m2', ax=plt.subplot(gs[10]))

    obs_vals = df_stan_data.c_aft_m2.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='c_aft_m2', ax=plt.subplot(gs[11]))


    obs_vals = df_syn_stan_data.c_pre1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='c_pre1', ax=plt.subplot(gs[12]))

    obs_vals = df_syn_stan_data.c_aft1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='c_aft1', ax=plt.subplot(gs[13]))
    '''
    fig.patch.set_visible(False)

    if fout != None:
        plt.savefig(fout, format='svg')
    plt.show()


def plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_syn1', ax=None, p_title=None, xlabel_size=8,
                               title_size=8):
    # plots the rank of a bunch of values vs. the posterior samples:
    # ie. for each of the 279 lrr_syn_1_pred, what is the rank in the posterior predictive distribution?
    # want a uniform distribution

    df_syn_lrr1_pred = df_fit[[c for c in df_fit.columns if c.startswith(var_plot + '_pred[')]]

    df_syn_by_it_pred = df_syn_lrr1_pred.T

    # subtract the oberved values
    df_syn_by_it_pred_diff = df_syn_by_it_pred.subtract(obs_vals, axis=0)
    df_syn_by_it_pred_diff['frac_below_0'] = df_syn_by_it_pred_diff.apply(lambda r: calc_p_val_from_threshold(r),
                                                                          axis=1)  # apply row-wise

    if ax == None:
        plt.figure()
        plt.hist(df_syn_by_it_pred_diff['frac_below_0'])
        plt.xlabel('fraction below observed')
        plt.title('histogram of where in distribution the vars fall {}'.format(var_plot))
        plt.show()
    else:
        ax.hist(df_syn_by_it_pred_diff['frac_below_0'], color='black')
        ax.set_xlabel('fraction simulated below observed', size=xlabel_size)
        ax.set_xlim([0, 1])
        if p_title == None:
            ax.set_title('{} obs rank in simulated'.format(var_plot), size=title_size)
        else:
            ax.set_title(p_title, size=title_size)
        ax.patch.set_visible(False)
    return df_syn_by_it_pred_diff


def plot_ppc(ax, df_fit, var_plot, obs_data, t_func, p_title, xlabel, title_size=8):
    '''
    Plot posterior predictive checks for a particular summary statistic of observed and
    posterior simulations.

    for example: var_plot = c_pre1:
    take the observed mean of the c_pre1, and plot the histogram of the mean of c_pre1 across all the simualtions.

    Parameters
    ----------
    ax : plt.axes subclass
        to plot into
    df_fit : stan df_fit object
    var_plot : str,
        such as 'c_pre'
    obs_data : dict
        dictionary fed into stan model for in sm.sampling() containing the observed data var_plot as key
    t_func : function(), needs to be able to accept axis argument
        ie. np.mean()
    p_title : str
        plot title
    xlabel : str
        xlabel for plot

    '''
    print(t_func, p_title, xlabel)

    rep_cols = [c for c in df_fit.columns if c.startswith(var_plot + '_rep')]
    if rep_cols == []:
        rep_cols = [c for c in df_fit.columns if c.startswith(var_plot + '_pred')]
    if rep_cols == []:
        print('Dataframe does not contain _pred or _rep columns to plot PPCs.', df_fit.columns)
        raise

    c_pre_rep = df_fit[rep_cols]
    c_pre_rep_t = t_func(np.ma.masked_invalid(c_pre_rep), axis=1)
    obs_data_vals = obs_data[var_plot]
    obs_t = t_func(np.ma.masked_invalid(obs_data_vals))

    ax.hist(np.ma.masked_invalid(c_pre_rep_t),
            bins=20,
            alpha=0.5,
            label='MCMC samples',
            color='black')
    ax.axvline(obs_t, label='observed data', color='red')
    ax.legend(loc='upper left', labelspacing=0.1, prop={'size': 5})
    ax.set_title(p_title, size=title_size)
    p_value = np.sum(i > obs_t for i in c_pre_rep_t) / len(c_pre_rep_t)
    return p_value


def compare_posterior_std_synon(df_fit, stan_data, var_to_rank, obs_vals, ax=None, p_title=None, xlabel_size=8,
                                title_size=8):
    ''' plot the rank of the standard deviation of particular variables.

    :param df_m_curr: df_fit
    :param stan_data: dictionary-form data to feed into stan model
    :param var_to_rank: ie. lrr_m1_pred[
    :param obs_vals: the observed value
    '''
    # calculate the observed synonymous standard deviation
    var_cols = [c for c in df_fit.columns if c.startswith(var_to_rank)]
    df_fit_var = df_fit[var_cols]
    df_cod_by_it_var = df_fit_var.T  # so now is i_aa by n_mcmc_samples df

    # add amino acid idx, groupby it, and calculate stdev
    df_cod_by_it_var['s'] = stan_data['s'].values
    df_cod_std_var = df_cod_by_it_var.groupby('s').std()

    # subtract the observed stdev
    df_cod_std_var_diff = df_cod_std_var.subtract(obs_vals, axis=0)

    # calculate fraction of stdeviations from posterior that are smaller than the observed standard deviations
    df_cod_std_var_diff['frac_below_0'] = df_cod_std_var_diff.apply(lambda r: calc_p_val_from_threshold(r),
                                                                    axis=1)  # apply row-wise

    # this means that the stdeviation observed is in the posterior samples is smaller than the observed stdev
    if ax == None:
        plt.figure()
        plt.hist(df_cod_std_var_diff['frac_below_0'])
        plt.xlabel('fraction below observed')
        plt.title('distribution of stdevs of synonymous mutants for {}'.format(var_to_rank))
        plt.show()

    else:
        ax.hist(df_cod_std_var_diff['frac_below_0'], color='black')
        ax.set_xlabel('fraction below observed', size=xlabel_size)
        if p_title == None:
            ax.set_title('stdeviations {}'.format(var_to_rank), size=title_size)
        else:
            ax.set_title(p_title, size=title_size)
        ax.set_xlim([0, 1])

    return df_cod_std_var_diff


def ppc_figures_df_fit(df_fit,
                       stan_data, df_stan_data, df_syn_stan_data):
    '''

    # get the distributions

    plots and returns calculated posterior predictive checks for figures

    1) Can you capture the synonymous mutant deviations
        hist: rank of synonymous standard deviations of log read ratio rep 1, lrr1 and lrr2 for mutants
    2) Is there a bias in the inferred synonymous wt lrrs?
        hist: rank of synnonymous wt toxin lrr1_syn and lrr2_syn
    3) Is there a bias in the inferred lrrs of codon mutants?
        hist: rank of lrr1 codon fitness in simulated codon fitnesses, and lrr2, for synonymous and mutants
    4) Is there a bias in the inferred amino acid lrrs?
        hist: lrr amino acid fitness in simulated lrr amino acid fitnesses

    5) mean + stdeviation of synonymous wt toxins lrr vs. simulated mean
    6) mean + stdeviation of synonymous mutants toxins lrr vs. simulated mean

    7) to check that the pre-counts are fit well:
        is the mean and stdev of the observed c_pre1 and c_pre_m1 captured? (comparing negbin vs. lognormal-poisson)



    :param df_fit: expected to have log read ratios 'lrr'

    '''

    # rename for compatibility
    df_syn_stan_data = df_syn_stan_data.rename(columns={'lrr1': 'lrr_syn1', 'lrr2': 'lrr_syn2'})
    df_stan_data = df_stan_data.rename(columns={'lrr1': 'lrr_m1', 'lrr2': 'lrr_m2'})

    # 1)
    print('--------------------------------------------------------------------------------------')
    print('------- 1) Can you capture the synonymous mutant deviations?')
    obs_mutant_stds = df_stan_data.groupby('g_idx').std()
    df_cod_std_lrr1_diff_full = compare_posterior_std_synon(df_fit, stan_data, 'lrr_m1_pred[',
                                                            obs_mutant_stds['lrr_m1'])
    df_cod_std_lrr2_diff_full = compare_posterior_std_synon(df_fit, stan_data, 'lrr_m2_pred[',
                                                            obs_mutant_stds['lrr_m2'])

    # 2) plot rank distribution of synonymous mutants
    print('--------- 2)Is there a bias in the inferred synonymous wt lrrs?')
    obs_vals = df_syn_stan_data.lrr_syn1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_syn1')

    obs_vals = df_syn_stan_data.lrr_syn2.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_syn2')

    # 3) plot rank distribution of mutant codons
    print('---------  3) Is there a bias in the inferred mutant codon lrrs?')
    obs_vals = df_stan_data.lrr_m1.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_m1')

    obs_vals = df_stan_data.lrr_m2.values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_m2')

    # 4)
    print('---------  4) Is there a bias in the inferred mutant amino acid mean lrrs?')
    obs_vals = stan_data['lrr_aa_m1'].values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_aa_m1')

    obs_vals = stan_data['lrr_aa_m2'].values
    plot_rank_distribution_var(df_fit, obs_vals, var_plot='lrr_aa_m2')

    # 5)
    print('---------  5) Is there an overall shift in the mean of synonymous wt log read ratios?')
    plot_all_ppcs(df_fit, 'lrr_syn1', df_syn_stan_data)
    plot_all_ppcs(df_fit, 'lrr_syn2', df_syn_stan_data)

    # 6)
    print('---------  6) Is there an overall shift in the mean of codon mutant log read ratios?')
    plot_all_ppcs(df_fit, 'lrr_m1', df_stan_data)
    plot_all_ppcs(df_fit, 'lrr_m2', df_stan_data)

    # 7)
    print('---------  7)  Can we recapitulate counts??')
    plot_all_ppcs(df_fit, 'c_pre1', df_syn_stan_data)
    plot_all_ppcs(df_fit, 'c_pre2', df_syn_stan_data)

    plot_all_ppcs(df_fit, 'c_pre_m1', df_stan_data)
    plot_all_ppcs(df_fit, 'c_pre_m2', df_stan_data)


def plot_all_ppcs(df_fit, var_plot, obs_data, plot_out=None):
    '''
    Plot posterior predictive distribution for mean, stdev, max and range of count data

    Parameters
    ----------
    df_fit : stan df_fit object
    var_plot : str,
        such as 'c_pre'
    obs_data : dict
        dictionary fed into stan model for in sm.sampling() containing the observed data var_plot as key
    '''

    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(2, 2)

    ax0 = plt.subplot(gs[0])
    p_val_mean = plot_ppc(ax0, df_fit, var_plot, obs_data, np.mean, 'mean',
                          'mean')
    ax1 = plt.subplot(gs[1])

    p_val_std = plot_ppc(ax1, df_fit, var_plot, obs_data, np.std, 'stdev ',
                         'stdev')

    ax2 = plt.subplot(gs[2])
    p_val_max = plot_ppc(ax2, df_fit, var_plot, obs_data, np.max, 'max ', 'max')

    ax3 = plt.subplot(gs[3])
    p_val_range = plot_ppc(ax3, df_fit, var_plot, obs_data, get_range, 'range ',
                           'range')

    if plot_out == None:
        plt.show()
    else:
        plt.savefig(plot_out, format='pdf')

    return p_val_mean, p_val_std, p_val_max, p_val_range


#################################################################################################
############## functions to validate against synthetic data
#################################################################################################


def get_negbin_data(data_all):
    # function to get negative binomial fits of the pre selection reads
    # for example to match the G62L high observed dataset in pre read counts
    nb_fit1 = statTools.get_negbin_fit(data_all['c_pre1'])
    nb_fit2 = statTools.get_negbin_fit(data_all['c_pre2'])

    nb1_n = nb_fit1[1]
    nb1_p = nb_fit1[2]
    nb2_n = nb_fit2[1]
    nb2_p = nb_fit2[2]

    # plot the fractional rank of the mean log read ratio stdeviation for synonymous mutants

    nb_fit_m1 = statTools.get_negbin_fit(data_all['c_pre_m1'])

    m1_n = nb_fit_m1[1]
    m1_p = nb_fit_m1[2]

    nb_fit_m2 = statTools.get_negbin_fit(data_all['c_pre_m2'])

    m2_n = nb_fit_m2[1]
    m2_p = nb_fit_m2[2]

    return (nb1_n, nb1_p, nb2_n, nb2_p), (m1_n, m1_p, m2_n, m2_p)


def correlate_inferred(df_perc, true_diff_gr, figsize=(2, 2)):
    # plot correlation of inferred amino acid growth rates, and true growth rates

    plt.figure(figsize=figsize)
    plt.scatter(true_diff_gr, df_perc.mean_fit)
    plt.plot([-1, 3], [-1, 3], c='black')
    plt.xlabel('true growth rate')
    plt.ylabel('inferred growth rate')
    plt.show()

    corr = st.pearsonr(true_diff_gr, df_perc.mean_fit)
    print('pearsonr', corr)
    return corr


def correlate_inferred_errorbar(df_perc, true_diff_gr, figsize=(2, 2), fout=None):
    # plot correlation of inferred amino acid growth rates, and true growth rates

    plt.figure(figsize=figsize)
    # plt.scatter(true_diff_gr, df_perc.mean_fit)
    plt.errorbar(true_diff_gr, df_perc.mean_fit,

                 # to make a (2,N) array with absolute deviations in either direciton
                 yerr=np.array([df_perc['2.5%'].values - df_perc['mean_fit'].values,
                                df_perc['mean_fit'].values - df_perc['97.5%'].values]),
                 elinewidth=0.3,
                 alpha=0.5,
                 # color='black',
                 fmt=".")
    offset = 0.4

    xmin = min(true_diff_gr) - offset
    xmax = max(true_diff_gr) + offset
    ymin = df_perc.mean_fit.min() - offset
    ymax = df_perc.mean_fit.max() + offset
    plt.ylim([ymin, ymax])
    plt.xlim([xmin, xmax])
    plt.plot([ymin, ymax], [ymin, ymax], c='black')

    plt.xlabel('true growth rate')
    plt.ylabel('inferred growth rate')
    if fout != None:
        plt.savefig(fout, format='svg')
    plt.show()

    corr = st.pearsonr(true_diff_gr, df_perc.mean_fit)
    print('pearsonr', corr)
    return corr


def corr_hdi_res(df_perc, true_diff_gr, figsize=(5, 5)):
    # correlate residual (true growth rate difference to inferred growth rate) vs. hdi length into the direction

    df_perc['resid'] = df_perc['mean_fit'] - true_diff_gr

    def get_correct_hdi_val(r):
        if r.resid > 0:
            return r['mean_fit'] - r['2.5%']
        else:
            return r['mean_fit'] - r['97.5%']

    df_perc['hdi_diff'] = df_perc.apply(lambda r: get_correct_hdi_val(r), axis=1)

    plt.figure(figsize=figsize)
    plt.scatter(df_perc.resid, df_perc.hdi_diff, s=1)
    plt.xlabel('inferred GR - true GR')
    plt.ylabel('HDI')
    plt.plot([-2, 1], [-2, 1], c='black')
    plt.show()

    corr = st.pearsonr(df_perc.resid, df_perc.hdi_diff)
    return corr


def calc_frac_in_hdi(df_perc, true_diff_gr, tail_prob=0.05):
    # calc how many of the true growth rates are within the percentile HDI

    # fetch the correct limits from df_percent
    upper_lim = '{:.1f}%'.format((1 - tail_prob / 2) * 100)
    lower_lim = '{:.1f}%'.format((tail_prob) / 2 * 100)

    # hacky to prevent 95.0% 5.0%
    if upper_lim[-2] == '0':
        upper_lim = '{:d}%'.format(int((1 - tail_prob / 2) * 100))
        lower_lim = '{:d}%'.format(int((tail_prob) / 2 * 100))

    # print(upper_lim, lower_lim)
    df_perc['in_hdi'] = (df_perc[upper_lim] > true_diff_gr) & (true_diff_gr > df_perc[lower_lim])

    fraction_hit = np.sum(df_perc.in_hdi) / len(df_perc)
    return fraction_hit


def calc_frac_in_hdi_cutoffs(df_perc, true_diff_gr):
    for p_tail in [0.05, 0.1, 0.2, 0.50]:
        in_hdi = calc_frac_in_hdi(df_perc, true_diff_gr, tail_prob=p_tail)
        print('{}% HDI has so many true values in the interval: {}'.format(1 - p_tail, in_hdi))


def get_synthetic_data(nb_params=[],
                       nb_params_mutants=[],
                       s=[],
                       true_gr=[],
                       true_gr_syn=None,
                       n_syn_codons=None,
                       read_ratios=1,
                       std_noise_reps=0,
                       std_noise_syn=0,
                       std_noise_syn_wt=0,
                       rand_seed=3):
    '''
    create synthetic data for 2 replicates.
    matched for
    :param nb_params: vector of negative binomial 1 n, nb1 p, nb2 n, nb2 p
    :param s: [1, 1, 2, 2, 2, ..., 55, 55]: vector of length codon number n, the values indicate which codon codes
    for which amino acid
    :param true_gr: list of true growth ratios, should be length n_aa
    :return:
    '''
    np.random.seed(rand_seed)

    if type(std_noise_syn) == list:
        print('using list of synonymous stdevs supplied')
    elif type(std_noise_syn) == float:
        std_noise_syn = std_noise_syn * np.ones(len(true_gr))
        print('broadcasting scalar synonymous noise')

    # generate a codon length list of true growth rates, with noise for synonymous mutants.
    true_gr_codons = [np.random.normal(true_gr[i - 1], std_noise_syn[i - 1]) for i in s]
    true_gr_codons_1 = [np.random.normal(gr, std_noise_reps) for gr in true_gr_codons]
    true_gr_codons_2 = [np.random.normal(gr + np.log(1 / read_ratios), std_noise_reps) for gr in true_gr_codons]

    true_gr_codons_syn = [np.random.normal(true_gr_syn, std_noise_syn_wt) for _ in range(n_syn_codons)]
    true_gr_codons_syn_1 = [np.random.normal(gr, std_noise_reps) for gr in true_gr_codons_syn]
    true_gr_codons_syn_2 = [np.random.normal(gr + np.log(1 / read_ratios), std_noise_reps) for gr in true_gr_codons_syn]

    # unpack nb params
    nb1_n, nb1_p, nb2_n, nb2_p = nb_params

    # if separate mutant negative binomial fits are given
    if nb_params_mutants != []:
        m1_n, m1_p, m2_n, m2_p = nb_params_mutants
    else:
        m1_n, m1_p, m2_n, m2_p = nb1_n, nb1_p, nb2_n, nb2_p

    n_codon_muts = len(s)
    n_aa_muts = len(true_gr)

    # creating synonymous wt mutants
    a1 = nb1_n
    beta_scale1 = 1 / nb1_p - 1
    a2 = nb2_n
    beta_scale2 = 1 / nb2_p - 1

    a1_m = m1_n
    a2_m = m2_n
    beta_scale_m1 = 1 / m1_p - 1
    beta_scale_m2 = 1 / m2_p - 1

    # sampling
    f_pre1 = st.gamma.rvs(a1, scale=beta_scale1, size=n_syn_codons)
    # f_pre1 = st.nbinom.rvs(nb1_n, nb1_p, size=n_syn_codons)
    c_pre1 = st.poisson.rvs(f_pre1)
    r_syn1 = np.exp(true_gr_codons_syn_1)
    f_aft1 = f_pre1 * r_syn1
    c_aft1 = st.poisson.rvs(f_aft1)

    f_pre2 = st.gamma.rvs(a2, scale=beta_scale2, size=n_syn_codons)
    # f_pre2 = st.nbinom.rvs(nb2_n, nb2_p, size=n_syn_codons)
    c_pre2 = st.poisson.rvs(f_pre2)
    r_syn2 = np.exp(true_gr_codons_syn_2)
    f_aft2 = f_pre2 * r_syn2
    c_aft2 = st.poisson.rvs(f_aft2)

    # creating nonsynonymous amino acid mutants
    n_codons = len(s)
    f_pre_m1 = st.gamma.rvs(a1_m, scale=beta_scale_m1, size=n_codons)
    # f_pre_m1 = st.nbinom.rvs(nb1_n, nb1_p, size=n_codons)
    c_pre_m1 = st.poisson.rvs(f_pre_m1)
    r = np.exp(true_gr_codons_1)
    f_aft_m1 = f_pre_m1 * r
    c_aft_m1 = st.poisson.rvs(f_aft_m1)

    f_pre_m2 = st.gamma.rvs(a2_m, scale=beta_scale_m2, size=n_codons)
    # f_pre_m2 = st.nbinom.rvs(nb2_n, nb2_p, size=n_codons)
    c_pre_m2 = st.poisson.rvs(f_pre_m2)
    r = np.exp(true_gr_codons_2)
    f_aft_m2 = f_pre_m2 * r
    c_aft_m2 = st.poisson.rvs(f_aft_m2)

    data_synthetic = {
        'K'       : int(len(c_pre1)),
        'c_pre1'  : c_pre1,
        'c_aft1'  : c_aft1,
        'c_pre2'  : c_pre2,
        'c_aft2'  : c_aft2,
        'Km'      : n_aa_muts,
        'N'       : n_codon_muts,
        's'       : s,
        'c_pre_m1': c_pre_m1,
        'c_aft_m1': c_aft_m1,
        'c_pre_m2': c_pre_m2,
        'c_aft_m2': c_aft_m2
    }
    return data_synthetic, true_gr, true_gr_syn
