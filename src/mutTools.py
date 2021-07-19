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

