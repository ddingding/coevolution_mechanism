# classify reads and call single mutants

import yaml
import pandas as pd
from pipelineTools import class_samples

# read config files
config_dics = yaml.safe_load(open("./ex51_config.yaml", "r"))
df_config = pd.read_csv(open("./ex51_config.csv", "r"))


if __name__ == "__main__":
    # for hiseq run 1
    # din = '/n/groups/marks/users/david/ex51/04_split_by_ats/'
    # dout = '/n/groups/marks/users/david/ex51/05_class/'

    # for hiseq run 2
    # din = '/n/groups/marks/users/david/ex51/04_split_by_ats_4021W/'
    # dout = '/n/groups/marks/users/david/ex51/05_class_4021W/'

    # for mcs hiseq run1
    # din = '/n/groups/marks/users/david/ex51/03filtered/mcs/'
    # dout = '/n/groups/marks/users/david/ex51/05_class/mcs/'

    # for mcs hiseq run 2
    din = "/n/groups/marks/users/david/ex51/03filtered_4021W/mcs/"
    dout = "/n/groups/marks/users/david/ex51/05_class_4021W/mcs/"

    class_samples(din, dout, df_config)
