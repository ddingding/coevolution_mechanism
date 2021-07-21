'''
This script creates the directory structure for processing of raw read files.


'''

import os
from constants import DL_DIR, T_SINGLE_1_DIRS, T_SINGLE_2_DIRS, AT_SINGLE_DIRS, AT_COMBO_DIRS, \
    AT_COMBO_RAW_DIR, \
    AT_SINGLE_RAW_DIR_M1, AT_SINGLE_RAW_DIR_M2, \
    T_SINGLE_1_RAW_DIR_M1, T_SINGLE_1_RAW_DIR_M2, \
    T_SINGLE_2_RAW_DIR_M1, T_SINGLE_2_RAW_DIR_M2
from exp_to_files import at_single_machine1_files, at_single_machine2_files, \
    at_combo_files, \
    t_muts_1_machine1_files, t_muts_1_machine2_files, \
    t_muts_2_machine1_files, t_muts_2_machine2_files

########################################################################################################################
# create directories for processing
for d in T_SINGLE_1_DIRS + T_SINGLE_2_DIRS + AT_SINGLE_DIRS + AT_COMBO_DIRS:
    os.system('mkdir {}'.format(d))

#########################################################################################################################
# moving raw fastq files into the respective directories



for f in at_combo_files:
    os.system('mv {}{} {}'.format(DL_DIR, f, AT_COMBO_RAW_DIR))

for f in at_single_machine1_files:
    os.system('mv {}{} {}'.format(DL_DIR, f, AT_SINGLE_RAW_DIR_M1))

for f in at_single_machine2_files:
    os.system('mv {}{} {}'.format(DL_DIR, f, AT_SINGLE_RAW_DIR_M2))

for f in t_muts_1_machine1_files:
    os.system('mv {}{} {}'.format(DL_DIR, f, T_SINGLE_1_RAW_DIR_M1))

for f in t_muts_1_machine2_files:
    os.system('mv {}{} {}'.format(DL_DIR, f, T_SINGLE_1_RAW_DIR_M2))

for f in t_muts_2_machine1_files:
    os.system('mv {}{} {}'.format(DL_DIR, f, T_SINGLE_2_RAW_DIR_M1))

for f in t_muts_1_machine2_files:
    os.system('mv {}{} {}'.format(DL_DIR, f, T_SINGLE_2_RAW_DIR_M2))