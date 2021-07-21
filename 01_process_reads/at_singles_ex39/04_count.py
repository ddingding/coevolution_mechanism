import pandas as pd
from mutTools import get_f_pairs_sample
from constants import AT_SINGLE_CALLED_DIR_BOTH, T_SINGLE_2_COUNTED_DIR

class_din = AT_SINGLE_CALLED_DIR_BOTH
count_out = AT_SINGLE_COUNTED_DIR
template = 'pard'

#df_config = pd.read_csv(open("./ex47_config.csv", "r"))
# get the configuration...

#write files out
sample_to_df_counts = get_f_pairs_sample(class_din, df_config, dout = count_out,template=template)