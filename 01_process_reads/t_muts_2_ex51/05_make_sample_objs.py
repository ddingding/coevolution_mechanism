# make sample objects and timecourse objects out of classified reads.
from os import listdir
import yaml
import pandas as pd


# read config files
config_dics = yaml.safe_load(open("./ex51_config.yaml", "r"))
# py2
# df_config = pd.read_csv(open('./ex51_config.csv', 'r'))

# py3
df_config = pd.read_csv(open("./ex51_config.csv", "r", encoding="ISO-8859-1"))

# class_din = '/n/groups/marks/users/david/ex51/05_class_concat/'
class_din = "/n/groups/marks/users/david/ex51/05_class_concat/mcs_renamed/"
pickle_out = "/n/groups/marks/users/david/ex51/06_pickles/"
plot_out = "/n/groups/marks/users/david/ex51/07_plots/"

#making samples objects. had to do this 4 times to get through all the sample opjects, memory issues.
dic_samples = make_all_sample_obj(class_din, pickle_out, df_config)


s_name_list = [
    f.rstrip("_sample_obj.p")
    for f in listdir(pickle_out)
    if f.endswith("_sample_obj.p")
]

# making all the tc objects and pickling
print("s_name_list", sorted(s_name_list))
tc_list = make_all_possible_tcs_single(s_name_list, pickle_out, df_config)
