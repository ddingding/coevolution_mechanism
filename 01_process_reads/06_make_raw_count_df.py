# this script converts a timecourse object to a dataframe containing read counts pre- and post-selection
# and this data is used afterwards for the bayesian inference.

from os import listdir
import pickle

import libClassTools as lct


def get_all_raw_count_dfs(pickle_out, df_out):
    tc_list = [pickle_out + f for f in listdir(pickle_out) if f.endswith('tc.p')]

    tot_num_tcs = len(tc_list)
    print('number of timecourse objects to convert to df: {}'.format(len(tc_list)))

    for tc_path in tc_list:
        print('trying ' + pickle_out + tc_path)
        tc = pickle.load(open(tc_path, 'rb'), encoding='latin1')
        fname = df_out + tc_path.split('/')[-1][:-2] + '.csv'
        lct.get_raw_counts_tc(tc, fname)
        print('finished ' + fname)


# for antitoxin single mutant library timecourses
pickle_out = '/Users/davidd/DropboxLinks/DropboxHMS/parESingleLibrary/ex39_libraryRun/illumina/pickles/'
df_out = '/Users/davidding/PycharmProjects/pareSingleLibrary2/codebase/pairedEnd/x51/df/tc/ex39/'
get_all_raw_count_dfs(pickle_out, df_out)

# for t_muts_1 timecourses
pickle_out = '/n/groups/marks/users/david/ex47/06_pickles/pickles_concat/'
df_out = '/n/groups/marks/users/david/ex47/06_pickles/df/'
get_all_raw_count_dfs(pickle_out, df_out)

# for t_muts_2 timecourses
pickle_out = '/Users/davidding/Dropbox (HMS)/parESingleLibrary/ex51_set_up_additional_mutants/illumina/pickles/'
df_out = '/Users/davidding/PycharmProjects/pareSingleLibrary2/codebase/pairedEnd/x51/df/tc/ex51/'
get_all_raw_count_dfs(pickle_out, df_out)

