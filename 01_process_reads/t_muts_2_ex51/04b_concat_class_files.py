# merge the classified read output from 2 sequencing runs
# takes 2 dirs, and concatenates all the files with same names
from pipelineTools import concat_files, concat_fs_in_dirs

path1 ='/n/groups/marks/users/david/ex51/05_class/'
path2='/n/groups/marks/users/david/ex51/05_class_4021W/'

dout ='/n/groups/marks/users/david/ex51/05_class_concat/'
concat_fs_in_dirs(path1, path2, dout)
