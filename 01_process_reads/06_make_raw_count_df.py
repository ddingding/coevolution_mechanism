from os import listdir
import libClassTools as lct
import pickle


'''
# for ex47 timecourses 
pickle_out = '/n/groups/marks/users/david/ex47/06_pickles/pickles_concat/'
df_out = '/n/groups/marks/users/david/ex47/06_pickles/df/'
'''

# for ex51 timecourses
pickle_out = '/Users/davidding/Dropbox (HMS)/parESingleLibrary/ex51_set_up_additional_mutants/illumina/pickles/'
df_out = '/Users/davidding/PycharmProjects/pareSingleLibrary2/codebase/pairedEnd/x51/df/tc/ex51/'

if __name__ == '__main__':
    tc_list = [pickle_out + f for f in listdir(pickle_out) if f.endswith('tc.p')]

    tot_num_tcs = len(tc_list)
    print('number of timecourse objects to convert to df: {}'.format(len(tc_list)))

    for tc_path in tc_list:
        print('trying ' + pickle_out + tc_path)
        tc = pickle.load(open(tc_path, 'rb'), encoding='latin1')
        fname = df_out + tc_path.split('/')[-1][:-2] + '.csv'
        lct.get_raw_counts_tc(tc, fname)
        print('finished ' + fname)
