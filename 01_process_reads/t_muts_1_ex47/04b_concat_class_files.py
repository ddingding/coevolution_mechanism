#concatenate machine 1 and machine 2 classified read output

from pipelineTools import concat_fs_in_dirs
from constants import  T_SINGLE_1_CALLED_DIR_M1, T_SINGLE_1_CALLED_DIR_M2, T_SINGLE_1_CALLED_DIR_BOTH
path1 =T_SINGLE_1_CALLED_DIR_M1
path2=T_SINGLE_1_CALLED_DIR_M2
dout =T_SINGLE_1_CALLED_DIR_BOTH

concat_fs_in_dirs(path1, path2, dout)
