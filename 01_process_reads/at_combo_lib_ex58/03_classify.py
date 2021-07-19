# take merged and filtered fasta files and call combinatorial mutants

from os import listdir
from os.path import isfile, join
from constants import AT_COMBO_FILTERED_DIR, AT_COMBO_CALLED_DIR
fasta_dir = AT_COMBO_FILTERED_DIR
output_dir = AT_COMBO_CALLED_DIR

fastas = [f for f in listdir(fasta_dir) if isfile(join(fasta_dir, f)) and
                      	f.endswith('fasta')]

for f in fastas:
	fin = fasta_dir + f
	fout = output_dir + f[:-6] + '.csv'

    my_cmd = 'process_one_combi.py {} {} '.format(fin, fout)
	print(my_cmd)
	os.system(my_cmd)

