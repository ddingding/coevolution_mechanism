'''
calls combinatorial antitoxin mutants for a fasta file that contains read information
'''

from mutTools import fasta_iter_py3, translate, hamming, sep_mutstr_into_codon_mutstr, mutstr_to_dic, mutdic_to_mutcodon
from constants import WT_PARD_DNA, WT_PARD_DNA_AA
import sys


fa_in = sys.argv[3]
file_out = sys.argv[4]

def process_combi_fasta_f(fa_f, file_out):
    with open(file_out, 'w') as fout:

        for rec in fasta_iter_py3(fa_f):
            n = rec.header
            read = rec.sequence

            # check that read is longer than a full length antitoxin
            if len(read) < 282:
                write_list = [n, 'len_read_less_282bp']
                fout.write('\t'.join(write_list) + '\n')
                continue

            # get the starting position of the antitoxin ORF from the sequence
            orf_start_pos = read.find('ATGGCAAACGTG')
            # if not found
            if orf_start_pos == -1:
                write_list = [n, 'noATGfound']
                fout.write('\t'.join(write_list) + '\n')
                continue


            orf_seq = read[orf_start_pos:]

            # check that orf is length 282
            if len(orf_seq) < 282:
                write_list = [n, 'len_orf_less_282bp']
                fout.write('\t'.join(write_list) + '\n')
                continue

            if len(orf_seq) > 282:
                write_list = [n, 'len_orf_more_282bp']
                fout.write('\t'.join(write_list) + '\n')
                continue

            # check that orf ends in stop codon
            if orf_seq[-3:] != 'TAG':
                write_list = [n, 'orf_end_not_stop']
                fout.write('\t'.join(write_list) + '\n')
                continue

            # now should only be full- length genes, starting with ATG, and ending in TAG

            # check for sequences that have way to many mutations (more than 9 nucleotide changes)
            if hamming(orf_seq, WT_PARD_DNA) > 9:
                write_list = [n, 'more_than_9_bp_muts']
                fout.write('\t'.join(write_list) + '\n')
                continue

            # catch all the wt antitoxin mutants
            if hamming(orf_seq, WT_PARD_DNA) == 0:
                write_list = [n, 'right', 'wtAT']
                fout.write('\t'.join(write_list) + '\n')
                continue

            # get mut str.
            muts = {}
            c = 0
            for (wtC, mutC) in zip(WT_PARD_DNA, orf_seq):
                if wtC != mutC:
                    muts[int(c)] = mutC
                c += 1

            codon_mut_aapos = set([int(int(i) / 3) for i in muts.keys()]) # set of mutations positions as aa positions
            num_codon_muts = len(codon_mut_aapos)

            # if number of codons mutated are more than 3
            if num_codon_muts > 3:
                #print(n, codon_mut_aapos)
                write_list = [n, 'more_than_3_codon_muts']
                fout.write('\t'.join(write_list) + '\n')
                continue

            # if codon mutations are not in the target positions D60, K63, E79.
            right_codons = [int(aa_pos_mut) in [60,63,79] for aa_pos_mut in codon_mut_aapos]
            if False in right_codons:
                write_list = [n, 'at_least_one_codon_mut_in_wrong_pos']
                fout.write('\t'.join(write_list) + '\n')
                continue


            mut_str = ",".join(str(n) + ":" + mutC for n, mutC in muts.items())

            # make list of mut_codons
            codon_muts = []
            aa_pos_to_muts = sep_mutstr_into_codon_mutstr(mut_str)
            for aa_pos, codon_mut_str in aa_pos_to_muts.items():
                wt_codon = WT_PARD_DNA[aa_pos * 3: aa_pos*3 + 3]
                codon_mut_dic = mutstr_to_dic(codon_mut_str)
                mut_codon = mutdic_to_mutcodon(codon_mut_dic, wt_codon)
                codon_muts.append(wt_codon + str(aa_pos * 3) + mut_codon)

            #make list of aa mutants
            aa_muts = []
            for codon_mut in codon_muts:
                #print(n, codon_mut)

                aa_pos = int(int(codon_mut[3:-3])/3)
                #print(n, aa_pos)
                wt_aa = WT_PARD_DNA_AA[aa_pos]
                mut_aa = translate(codon_mut[-3:])
                aa_muts.append(wt_aa + str(aa_pos) + mut_aa)
            print(codon_mut_aapos)
            write_list = [n, 'right', ':'.join(codon_muts), ':'.join(aa_muts)]
            fout.write('\t'.join(write_list) + '\n')
            continue

# for testing locally on laptop
#process_combo_fasta_f('./test.fa', open('./test_fa_out.text', 'w'))

process_combi_fasta_f(fa_in, file_out)
