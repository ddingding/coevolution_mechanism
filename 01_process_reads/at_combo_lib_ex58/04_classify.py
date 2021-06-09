
from os import listdir
from os.path import isfile, join
import pandas as pd
import pickle as p
indir = '/n/groups/marks/users/david/ex58/03_called/'
outdir = '/n/groups/marks/users/david/ex58/04_counted/'



fastas = [f for f in listdir(indir) if isfile(join(indir, f)) and
                      	f.endswith('.csv')]

'''
indir= './'
fastas = ['tmp.txt']
outdir = './'
'''

for f in fastas:
	print('doing file {}'.format(f))
	aa_count_dic = {}
	codon_count_dic = {}
	with open(indir +f, 'r') as fin:
		for l in fin:
			line_list = l.strip().split('\t')
			if line_list[1] == 'right':
				if line_list[2] == 'wtAT':
					if 'wtAT' in aa_count_dic:
						aa_count_dic['wtAT'] +=1
						codon_count_dic['wtAT'] += 1

					else:
						aa_count_dic['wtAT'] = 1
						codon_count_dic['wtAT'] = 1

				aa_muts =line_list[-1]
				if aa_muts in aa_count_dic:
					aa_count_dic[aa_muts] += 1
				else:
					aa_count_dic[aa_muts] = 1

				codon_muts =line_list[-2]
				if codon_muts in codon_count_dic:
					codon_count_dic[codon_muts] += 1
				else:
					codon_count_dic[codon_muts] = 1

	print(aa_count_dic)
	df_aa = pd.DataFrame([aa_count_dic])
	df_aa.to_csv(outdir + f[:-4] + '_aa_counts.csv')
	print(codon_count_dic)
	df_codon = pd.DataFrame([codon_count_dic])
	df_codon.to_csv(outdir + f[:-4] + '_codon_counts.csv')

	#with open(outdir + f[:-4] + '_aa_counts.csv', 'w') as f:
	#	print(aa_count_dic, file=f)

	#with open(outdir + f[:-4] + '_codon_counts.csv', 'w') as f:
	#	print(codon_count_dic, file=f)

	p.dump(aa_count_dic, open(outdir + f[:-4] + '_aa_counts.p', 'wb'))
	p.dump(codon_count_dic, open(outdir + f[:-4] + '_codon_counts.p', 'wb'))