# functions to manipulate libclass objects
import numpy as np
from copy import deepcopy
from math import sqrt, pow, exp
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns


def get_raw_counts_tc(tc, fout):
    '''
    to make the count dfs for a  timecourse object.
    fout is the name of the output file
    '''
    s_pre = tc.samples[0]
    s_post = tc.samples[1]
    df1 = make_count_df(s_pre)
    df2 = make_count_df(s_post)
    df_merge = merge_codon_df({"pre": df1, "post": df2})
    df_merge.to_csv(fout)


def make_count_df(s_obj):
    '''
    create a dataframe of read counts for each single mutant in a sample object
    '''
    cols = [
        "wt_aa",
        "aa_pos",
        "mut_aa",
        "wt_codon",
        "codon_pos",
        "mut_codon",
        "fit_codon",
        "fit_aa",
        "raw_count",
    ]

    df_codon = pd.DataFrame(columns=cols)

    for mkey, mut_obj in s_obj.sm_objects.items():
        fit_codon = mut_obj.fitness
        fit_aa = mut_obj.fitness_aa
        raw_count = mut_obj.rawCount

        try:
            scm_obj = list(mut_obj.codon_muts.values())[0]
        except AttributeError:
            scm_obj = list(mut_obj.codonMutants.values())[0]

        wt_aa = scm_obj.wtAA
        aa_pos = scm_obj.AA_pos
        mut_aa = scm_obj.mutAA

        wt_codon = scm_obj.wtCodon
        codon_pos = scm_obj.codonPos
        mut_codon = scm_obj.mutCodon

        row = pd.Series(
            [
                wt_aa,
                aa_pos,
                mut_aa,
                wt_codon,
                codon_pos,
                mut_codon,
                fit_codon,
                fit_aa,
                raw_count,
            ],
            cols,
        )
        df_codon = df_codon.append(row, ignore_index=True)

    return df_codon

def merge_codon_df(m_to_df):
    # take a dictionary of AT_mut to codon_df, like {'G62Y': df_codon, 'G62A':df_codon2}
    # and merge
    cols_to_merge_on = [
        "wt_aa",
        "aa_pos",
        "mut_aa",
        "wt_codon",
        "codon_pos",
        "mut_codon",
    ]
    df_merge = pd.DataFrame()
    for m, df in m_to_df.items():
        if df_merge.empty:
            mutkey1 = m
            df_merge = df
        else:
            df_merge = df_merge.merge(
                df,
                how="outer",
                left_on=cols_to_merge_on,
                right_on=cols_to_merge_on,
                suffixes=("", "_" + m),
            )

    # rename the columns 'fit' which came from the first df to 'fit_mutkey of first df'
    df_merge = df_merge.rename(
        index=str,
        columns={"fit_codon": "fit_codon_" + mutkey1, "fit_aa": "fit_aa_" + mutkey1},
    )
    return df_merge


# add mutkey column
def add_mutkey(df):
    df["mutkey"] = (
            df["wt_aa"].astype(str)
            + df["aa_pos"].astype(int).astype(str)
            + df["mut_aa"].astype("str")
    )
    return df


def add_codonkey(df):
    df['codonkey'] = (
            df["wt_codon"].astype(str)
            + df["codon_pos"].astype(int).astype(str)
            + df["mut_codon"].astype("str")
    )
    return df


def merge_df_counts(pin, pin2):
    '''
    expects path to 2 count files, creates mutkey and codonkey column, merges them based on codonkey.
    used for merging replicate codon mutant count dataframes.
    '''
    df = pd.read_csv(pin)
    df = add_mutkey(df)
    df = add_codonkey(df)
    # df_syn = df.loc[df["mutkey"].str[0] == df["mutkey"].str[-1]]

    df2 = pd.read_csv(pin2)
    df2 = add_codonkey(df2)
    # df2_syn = df.loc[df["mutkey"].str[0] == df["mutkey"].str[-1]]
    df_reps = df.merge(df2, left_on='codonkey', right_on='codonkey', suffixes=('_r1', '_r2'))
    return df_reps



