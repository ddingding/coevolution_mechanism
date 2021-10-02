import pandas as pd
from mutTools import get_f_pairs_sample
from constants import T_SINGLE_2_CALLED_DIR_NO_AT_BOTH, T_SINGLE_2_COUNTED_DIR

class_din = T_SINGLE_2_CALLED_DIR_NO_AT_BOTH
count_out = T_SINGLE_2_COUNTED_DIR
template = 'pare'

df_config = pd.read_csv(open("./ex51_config.csv", "r"))
df_config = df_config[:28]

#write files out
sample_to_df_counts = get_f_pairs_sample(class_din, df_config, dout = count_out,template=template)
