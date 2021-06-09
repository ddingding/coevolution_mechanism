#takes 2 dirs, and concatenates all the files with same names
from pipelineTools import concat_files, concat_fs_in_dirs
import os
from os import listdir
from os.path import isfile, join


path1 ='/n/groups/marks/users/david/ex47/05_class/'
path2='/n/groups/marks/users/david/ex47/05_class_3786W/'

dout ='/n/groups/marks/users/david/ex47/05_class_concat/'

if __name__ == '__main__':

    #rename accidentally putting files in wrong place
    din ='/n/groups/marks/users/david/ex47/'
    fs = [f for f in listdir(din) if isfile(join(din, f)) and
                f.endswith('.tsv')]
    for f in fs:
        old_fname = f
        new_fname = f.lstrip('05_class_concat')
        cmd = 'mv {0}{1} {2}{3}'.format(din, old_fname, din, new_fname)
        print(cmd)
        os.system(cmd)
    #concat_fs_in_dirs(path1, path2, dout)
