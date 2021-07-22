# mutTools
# set of functions to use for manipulating biological sequences, fasta files, mutation strings, ...

import numpy as np
import operator
import matplotlib.pyplot as plt
from itertools import groupby

from constants import CODON_TABLE


def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def complement(inbase):
    cDic = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return cDic[inbase]

def rev_complement(instr):
    compl = [complement(c) for c in instr]
    return "".join(compl[::-1])

def translate(inCodon):
    return CODON_TABLE[inCodon]

def fasta_iter(fasta_name):
    """ works for python 2 only
    given a fasta file,  yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def fasta_iter_py3(fasta_name):
    """
	given a fasta file,  yield tuples of header, sequence
	https://github.com/GuyAllard/fasta_iterator/blob/master/fasta_iterator/__init__.py
    """
    rec = None
    for line in open(fasta_name, "r"):
        if line[0] == ">":
            if rec:
                yield rec
            rec = FastaRecord(line.strip()[1:])
        else:
            rec.sequence += line.strip()

    if rec:
        yield rec


class FastaRecord:
    """
    Fasta Record, contains a header and a sequence
	# from https://github.com/GuyAllard/fasta_iterator/blob/master/fasta_iterator/__init__.py
    """

    def __init__(self, header="", sequence=""):
        self.header = header
        self.sequence = sequence

    def __str__(self):
        return ">{}\n{}".format(self.header, self.sequence)

def mutstr_to_dic(mutStr):
    # take mutstr like '4:T,5:A,8:A' and make a dictionary of position to mutant nt
    return dict(
        zip(
            [int(x.split(":")[0]) for x in mutStr.split(",")],
            [x.split(":")[1] for x in mutStr.split(",")],
        )
    )

def mutdic_to_mutcodon(mutpos_to_base, wtCodon):
    # take  {'100': 'A', '101': 'G'} and convert to mutcodon, which is a 3 letter base string.
    mutCodon = ""
    codonPos = min(mutpos_to_base.keys()) // 3 * 3

    for i in range(3):
        if i in [int(x) % 3 for x in mutpos_to_base]:
            mutCodon += mutpos_to_base[codonPos + i]
        else:
            mutCodon += wtCodon[i]
    return mutCodon


def mutStrToTupleList(mutStr):
    """
	gets a mutstring like 100:C,101:A
	and returns a sorted list of tuples [(100,'C'), (101, 'A')]
	"""
    mut_list = mutStr.split(",")
    return sorted(
        [(int(mut.split(":")[0]), mut.split(":")[1]) for mut in mut_list],
        key=operator.itemgetter(0),
    )

def sep_mutstr_into_codon_mutstr(full_mut_str):
    """
	get a mutstr like '4:T,5:A,8:A' which indexes a codon position and the base mutation
	and return a list of single codon mutStr
	ie. ['4:T,5:A', '8:A'] for 2 codon mutants.

	or if you get '4:T,5:A'
	should only return a single codon ['4:T,5:A']
	"""

    pos_to_mut_dic = dict(mutStrToTupleList(full_mut_str))  # also gets sorted
    aaPosToMuts = {}
    for p in pos_to_mut_dic.keys():
        aaPos = int(int(p) / 3)
        if aaPos in aaPosToMuts:
            aaPosToMuts[aaPos] = (
                aaPosToMuts[aaPos] + "," + str(p) + ":" + str(pos_to_mut_dic[p])
            )
        else:
            aaPosToMuts[aaPos] = str(p) + ":" + str(pos_to_mut_dic[p])
    return aaPosToMuts


# read the single mutant count strings into dictionaries
def count_class_file(class_fin, template, read_limit=1000000000):
    # read all the individual lines and reads with their mutation string descriptions
    c = 0
    mut_type_count = defaultdict(int)
    mut_count = defaultdict(int)
    wt_count = defaultdict(int)

    with open(class_fin, "r") as fin:
        csvIn = csv.reader(fin, delimiter="\t")
        for l in csvIn:
            mut_type = l[1]

            # count the number of different classified mutations
            mut_type_count[mut_type] += 1

            # deal with the codons mutants
            if mut_type.endswith("_codons"):
                strand = l[2]
                mut_str = l[5]
                mut_count[mut_str] += 1

            # count the wt dna
            if mut_type.startswith("wt_dna_"):
                wt_count["wt"] += 1

            c += 1
            if c > read_limit:
                break

    total_reads = sum(mut_type_count.values())

    # create a dataframe from this dictionary
    cols = [
        "wt_aa",
        "aa_pos",
        "mut_aa",
        "wt_codon",
        "codon_pos",
        "mut_codon",
        "raw_count"
    ]

    df_codon = pd.DataFrame(columns=cols)

    # convert the single mutant counts into codon mutants
    for mut_str, count in mut_count.items():
        dic_codon_mutstr = sep_mutstr_into_codon_mutstr(mut_str)
        if len(dic_codon_mutstr) == 1:

            # dict to store position to the mutated base
            mutpos_to_base = mutstr_to_dic(mut_str)

            # double check it's indeed a single mutant
            mutpos_codon_list = [int(i) // 3 for i in mutpos_to_base]
            mutpos_codon_list = mutpos_codon_list
            assert (
                    len(set(mutpos_codon_list)) == 1
            ), "mutation from mut_str found not in the same codon" + str(mut_str)

            AA_pos = int(mutpos_codon_list[0])

            codon_pos = AA_pos * 3

            if template == "pare":
                wt_codon = WT_PARE_DNA[codon_pos: codon_pos + 3]
            elif template == "pard":
                wt_codon = WT_PARD_DNA[codon_pos: codon_pos + 3]

            mut_codon = mutdic_to_mutcodon(mutpos_to_base, wt_codon)

            wt_aa = translate(wt_codon)
            mut_aa = translate(mut_codon)

            row = pd.Series(
                [
                    wt_aa,
                    AA_pos,
                    mut_aa,
                    wt_codon,
                    codon_pos,
                    mut_codon,
                    count,
                ],
                cols,
            )
            df_codon = df_codon.append(row, ignore_index=True)
    df_codon['codon_mutant'] = df_codon.wt_codon + df_codon.codon_pos.astype(int).astype(str) + df_codon.mut_codon
    df_codon['aa_mutant'] = df_codon.wt_aa + df_codon.aa_pos.astype(int).astype(str) + df_codon.mut_aa
    return df_codon


def get_counts_sample(class_fin, class_fin_post, template='pare'):
    #make a single dataframe that contains read counts pre and post selection
    df_pre = count_class_file(class_fin, template)
    df_post = count_class_file(class_fin_post, template)
    df_sample = df_pre.merge(df_post, left_on=['wt_codon', 'codon_pos', 'mut_codon'],
                             right_on=['wt_codon', 'codon_pos', 'mut_codon'], suffixes=('', '_post'))
    df_sample = df_sample.sort_values(by=['codon_pos', 'wt_aa', 'mut_aa'])
    df_sample = df_sample[['wt_aa', 'aa_pos', 'mut_aa', 'wt_codon', 'codon_pos', 'mut_codon',
                           'codon_mutant', 'aa_mutant', 'raw_count', 'raw_count_post']]
    return df_sample


def map_primer_to_gene_sample(primer_str, df_config):
    # for a primer number, find the gene and the sample number that it belongs to,
    # the sample number identifies the expression conditions
    g = df_config.loc[df_config.primer == int(primer_str), "gene"].iloc[0]
    sample = df_config.loc[df_config.primer == int(primer_str), "sample"].iloc[0]
    return str(g) + "_" + str(sample)


def get_f_pairs_sample(class_din, df_config, dout='./', template='pare', max_samples_read=100000000):
    # get a dictionary of pairs of files that belong to the same sample
    fs = [f for f in listdir(class_din) if f.endswith('class.tsv')]
    dic_id_to_fs = {}
    for f in fs:
        primer, at_mut, _ = f.split('_')
        gene_sample = map_primer_to_gene_sample(primer, df_config)
        sample_id = gene_sample + '_' + at_mut
        if sample_id in dic_id_to_fs:
            dic_id_to_fs[sample_id].append(f)
        else:
            dic_id_to_fs[sample_id] = [f]

    # create the respective count datafiles
    sample_to_df_counts = {}
    done_primer_at = []
    c = 0
    for sample_id in sorted(dic_id_to_fs.keys()):
        f_pairs = dic_id_to_fs[sample_id]
        # find the right files for timepoint t=0 and t=600
        t_to_f = dict(
            [(int(df_config.loc[df_config.primer.astype(int).astype(str) == x.split('_')[0]].t), x) for x in f_pairs])
        print(t_to_f)

        f1_name = t_to_f[0] # get the file that has the first timepoint
        f2_name = t_to_f[600]
        primer_at = f1_name[:-len('_class.tsv')]
        if primer_at not in done_primer_at:
            rep_suffix = '_rep1'
        else:
            rep_suffix = '_rep2'
        f_write_name = primer_at + rep_suffix + '.csv'
        if not isfile(dout + f_write_name):
            df_sample_counts = get_counts_sample(class_din +f1_name,
                                                 class_din +f2_name, template=template)
            sample_to_df_counts[primer_at + rep_suffix] = df_sample_counts
            df_sample_counts.to_csv(dout + f_write_name)
        c += 1
        if c == max_samples_read:
            break
    return sample_to_df_counts


# add mutkey column
def add_mutkey(df):
    df["mutkey"] = (
            df["wt_aa"].astype(str)
            + df["aa_pos"].astype(int).astype(str)
            + df["mut_aa"].astype("str")
    )
    return df


def add_codonkey(df):
    df['codonkey'] = (
            df["wt_codon"].astype(str)
            + df["codon_pos"].astype(int).astype(str)
            + df["mut_codon"].astype("str")
    )
    return df


def merge_df_counts(pin, pin2):
    '''
    expects path to 2 count files, creates mutkey and codonkey column, merges them based on codonkey.
    used for merging replicate codon mutant count dataframes.
    '''
    df = pd.read_csv(pin)
    df = add_mutkey(df)
    df = add_codonkey(df)
    # df_syn = df.loc[df["mutkey"].str[0] == df["mutkey"].str[-1]]

    df2 = pd.read_csv(pin2)
    df2 = add_codonkey(df2)
    # df2_syn = df.loc[df["mutkey"].str[0] == df["mutkey"].str[-1]]
    df_reps = df.merge(df2, left_on='codonkey', right_on='codonkey', suffixes=('_r1', '_r2'))
    return df_reps




