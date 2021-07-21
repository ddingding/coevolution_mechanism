import pandas as pd
from mutTools import get_f_pairs_sample
from constants import T_SINGLE_1_CALLED_DIR_BOTH, T_SINGLE_1_COUNTED_DIR

class_din = T_SINGLE_1_CALLED_DIR_BOTH
count_out = T_SINGLE_1_COUNTED_DIR
template = 'pare'

df_config = pd.read_csv(open("./ex47_config.csv", "r"))

#write files out
sample_to_df_counts = get_f_pairs_sample(class_din, df_config, dout = count_out,template=template)