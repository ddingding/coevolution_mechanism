
list_fs= []
for i in range(169, 217):
    f1 = 'D21-{}-4848L/210105Lau_D21-{}_1_sequence.fastq'.format(i, i)
    f2 = 'D21-{}-4848L/210105Lau_D21-{}_2_sequence.fastq'.format(i, i)

    list_fs.append(f1)
    list_fs.append(f2)

concat_fs_cmd = '\{' + ','.join(list_fs) + '\}'
pref = 'rsync -avP laub_ill@bmc-150.mit.edu:/mnt/bmc-pub15/Laub/210105Lau/'
final_cmd = pref + concat_fs_cmd + ' .'
print(final_cmd)


