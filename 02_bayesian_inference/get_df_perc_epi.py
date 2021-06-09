# this script creates a summary file for a MCMC run, including:
# what the different percentile highest density intervals for each double mutant are,
# what fraction of posterior double mutant mass overlaps the single mutant,
# and what fraction of the double mutant posterior overlaps the independent expectation value

# takes an MCMC output file, like df_sm_bhs_diff_w_aa_m_180_A16K.csv
# and a df that contains the independent expectations, like df_muts.csv
# and outputs a df_percentiles_epi csv file.

import sys
import stanTools
import pandas as pd

f = sys.argv[1] # the MCMC output file
pin = sys.argv[2] # the directory of the output files
fin_df_epi = sys.argv[3] # the .csv file containing the independent expectations for each double mutant.

sample_str = f[22:-4]  # like 'df_sm_bhs_diff_w_aa_m_180_A16K.csv'[22:-4]
fout_name = 'df_percentiles_epi' + sample_str + '.csv'

df_epi = pd.read_csv(fin_df_epi)
stanTools.create_df_perc_epi_p_log2(f, pin, fout_name, sample_str, df_epi, col_pre='diff_w_aa_m')
