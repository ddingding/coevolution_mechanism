# take merged, and filtered fasta files
# call mutants across entire length of gene
# write out per read

from os import listdir
from os.path import isfile, join

fasta_dir = '/n/groups/marks/users/david/ex58/02_filtered/'
output_dir = '/n/groups/marks/users/david/ex58/03_called/'


fastas = [f for f in listdir(fasta_dir) if isfile(join(fasta_dir, f)) and
                      	f.endswith('fasta')]

for f in fastas:

	fin = fasta_dir + f
	fout = output_dir + f[:-6] + '.csv'


    my_cmd = 'sbatch submitCombi.sh {} {} '.format(fin, fout)
	print(my_cmd)
	os.system(my_cmd)

