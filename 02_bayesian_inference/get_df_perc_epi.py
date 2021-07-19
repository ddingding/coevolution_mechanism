'''
this script creates a summary file for a bayesian output file, including:
what the different percentile highest density intervals for each double mutant are,
what fraction of posterior double mutant mass overlaps the single mutant,

and if the the nonlinear model has been run, and expected double mutant growth rates estimated, what fraction of the
double mutant posterior overlaps this independent expected value.

takes an MCMC output file, like df_sm_bhs_diff_w_aa_m_180_A16K.csv and summarizes for each mutant backgroud, the posterior distribtion in terms of percentiles.
and outputs a summary file

Optional step (commented out):
if the nonlinear model has been run, and double mutant expectation calculated, you can rerurn this script to also estimate how significantly different the double mutants grow from these expected growth rates
and a df that contains the independent expectations, like df_muts.csv, containing double mutant expectations,
#####################
'''
import sys
import stanTools
import pandas as pd

###########################################################
# create a summary file for the highest density intervals, and probability of growing differently from the antitoxin
# single mutant, for a MCMC output file
f = sys.argv[1]  # the MCMC output file
pin = sys.argv[2]  # the directory of the output files
sample_str = f[22:-4]  # f is like 'df_sm_bhs_diff_w_aa_m_180_A16K.csv'[22:-4]
fout_name = 'df_percentiles_diff_w_aa_m_' + sample_str + '.csv'
stanTools.create_df_perc(f, pin, fout_name, sample_str)

'''
############################################################
# optional step:
# to create a summary file that includes a column detailing significance of growing differently than the double
# mutant expectation
# this can only be run after 03_nonlinear_model has been run, and df_muts.csv containing the double mutant
# expectations has been created.
f = sys.argv[1]  # the MCMC output file
pin = sys.argv[2]  # the directory of the output files
fin_df_epi = sys.argv[3]  # the .csv file containing the independent expectations for each double mutant, ie. df_muts

sample_str = f[22:-4]
fout_name = 'df_percentiles_epi' + sample_str + '.csv'

df_epi = pd.read_csv(fin_df_epi)
stanTools.create_df_perc_epi_p_log2(f, pin, fout_name, sample_str, df_epi, col_pre='diff_w_aa_m')
'''