import pandas as pd
from mutTools import get_counts_sample
from constants import AT_SINGLE_CALLED_DIR_BOTH, AT_SINGLE_COUNTED_DIR

class_din = AT_SINGLE_CALLED_DIR_BOTH
count_out = AT_SINGLE_COUNTED_DIR
template = "pard"

f_at_rep1_t0 = "180716Lau_D18-6083_class.csv"
f_at_rep2_t0 = "180716Lau_D18-6084_class.csv"
f_at_rep1_t600 = "180716Lau_D18-6093_class.csv"
f_at_rep2_t600 = "180716Lau_D18-6094_class.csv"

print("counting replicate 1")
df_at_rep1 = get_counts_sample(
    class_din + f_at_rep1_t0, class_din + f_at_rep1_t600, template="pard"
)
df_at_rep1.to_csv(count_out + "at_singles_rep1.csv")

print("counting replicate 2")
df_at_rep2 = get_counts_sample(
    class_din + f_at_rep2_t0, class_din + f_at_rep2_t600, template="pard"
)
df_at_rep2.to_csv(count_out + "at_singles_rep2.csv")
