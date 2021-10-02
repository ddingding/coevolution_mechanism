# merge the classified read output from 2 sequencing runs
# takes 2 dirs, and concatenates all the files with same names
from pipelineTools import concat_files, concat_fs_in_dirs
from constants import T_SINGLE_2_CALLED_DIR_NO_AT_M1, T_SINGLE_2_CALLED_DIR_NO_AT_M2, T_SINGLE_2_CALLED_DIR_NO_AT_BOTH

path1 =T_SINGLE_2_CALLED_DIR_NO_AT_M1
path2= T_SINGLE_2_CALLED_DIR_NO_AT_M2

dout = T_SINGLE_2_CALLED_DIR_NO_AT_BOTH
concat_fs_in_dirs(path1, path2, dout)
