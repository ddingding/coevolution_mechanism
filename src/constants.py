
#######################################################################################################################
# Paths to set on your local machine
DL_DIR = './' # set this path to where the sra files have been downloaded to
WORK_DIR = './' # set this path to where you want the processed files to be

########################################################################################################################
# directories to be created by setup_dirs.py

# antitoxin combinatorial library
AT_COMBO_DIR = WORK_DIR + 'at_combo/' # master directory for antitoxin combinatorial library
AT_COMBO_RAW_DIR = AT_COMBO_DIR + '00_raw_fastq/' # raw read files
AT_COMBO_MERGED_DIR = AT_COMBO_DIR + '01_merged/' # for merged paired end reads
AT_COMBO_FILTERED_DIR = AT_COMBO_DIR + '02_filtered/' # for quality filtered reads
AT_COMBO_CALLED_DIR = AT_COMBO_DIR +  '03_called/' # for called reads, ie. which mutations are present in each read
AT_COMBO_COUNTED_DIR = AT_COMBO_DIR +  '04_counted/' # for counted files that summarize the called reads

# antitoxin single mutant library
AT_SINGLE_DIR = WORK_DIR + 'at_single/'
AT_SINGLE_RAW_DIR_M1 = AT_SINGLE_DIR + '01_raw_fastq_3272/'
AT_SINGLE_RAW_DIR_M2 = AT_SINGLE_DIR + '01_raw_fastq_3252/'
AT_SINGLE_MERGED_DIR_M1 =AT_SINGLE_DIR+ '02_merged_3272/'
AT_SINGLE_MERGED_DIR_M2 =AT_SINGLE_DIR+ '02_merged_3252/'
AT_SINGLE_FILTERED_DIR_M1 =AT_SINGLE_DIR+ '03_filtered_m1/'
AT_SINGLE_FILTERED_DIR_M2 =AT_SINGLE_DIR+ '03_filtered_m2/'
AT_SINGLE_CALLED_DIR_M1 = AT_SINGLE_DIR+ '04_called_mi1/'
AT_SINGLE_CALLED_DIR_M2 = AT_SINGLE_DIR+ '04_called_mi2/'
AT_SINGLE_CALLED_DIR_BOTH = AT_SINGLE_DIR+ '04_called_both/'
AT_SINGLE_COUNTED_DIR = AT_SINGLE_DIR+ '04_counted/'

# toxin single mutants in various antitoxin mutant backgrounds, batch1
T_SINGLE_1_DIR = WORK_DIR + 't_single_1/'

T_SINGLE_1_RAW_DIR_M1 = T_SINGLE_1_DIR + '00_fastq_raw_m1/'
T_SINGLE_1_RAW_DIR_M2 = T_SINGLE_1_DIR + '00_fastq_raw_3786_m2/'
T_SINGLE_1_CONCAT_DIR_M1 = T_SINGLE_1_DIR + '01_concat_m1/' # for concatenating read files from 2 lanes
T_SINGLE_1_CONCAT_DIR_M2 = T_SINGLE_1_DIR + '01_concat_m2/'
T_SINGLE_1_MERGED_DIR_M1 = T_SINGLE_1_DIR + '02_merged_m1/'
T_SINGLE_1_MERGED_DIR_M2 = T_SINGLE_1_DIR + '02_merged_m2/'
T_SINGLE_1_FILTERED_DIR_M1 = T_SINGLE_1_DIR + '03_filtered_m1/'
T_SINGLE_1_FILTERED_DIR_M2 = T_SINGLE_1_DIR + '03_filtered_m2/'
T_SINGLE_1_AT_SPLIT_DIR_M1 = T_SINGLE_1_DIR + '04_at_split_m1/'
T_SINGLE_1_AT_SPLIT_DIR_M2 = T_SINGLE_1_DIR + '04_at_split_m2/'
T_SINGLE_1_CALLED_DIR_M1 = T_SINGLE_1_DIR + '04_called_m1/'
T_SINGLE_1_CALLED_DIR_M2 = T_SINGLE_1_DIR + '04_called_m2/'
T_SINGLE_1_CALLED_DIR_BOTH = T_SINGLE_1_DIR + '04_called_both/'
T_SINGLE_1_COUNTED_DIR = T_SINGLE_1_DIR + '05_counted/'





# toxin single mutants in various antitoxin mutant backgrounds, batch2
T_SINGLE_2_DIR = WORK_DIR + 't_single_2/'

#######################################################################################################################
# to go from 3 letter amino acid code to one letter amino acid code
AA3_TO_AA1 = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}

# * signifies the initial barcoded wt toxin.
AA_LIST_ALPHABETICAL = "ACDEFGHIKLMNPQRSTVWY_*"

# standard codon table
# modified the stops to be '_'
# try ordering of CODONS by aa feature similarity
CODON_LIST_AA_CHEMISTRY = [
    ("TAA", "_"),
    ("TAG", "_"),
    ("TGA", "_"),
    ("GCG", "A"),
    ("GCA", "A"),
    ("GCC", "A"),
    ("GCT", "A"),
    ("ATA", "I"),
    ("ATC", "I"),
    ("ATT", "I"),
    ("CTG", "L"),
    ("CTA", "L"),
    ("CTC", "L"),
    ("CTT", "L"),
    ("TTG", "L"),
    ("TTA", "L"),
    ("ATG", "M"),
    ("GTG", "V"),
    ("GTA", "V"),
    ("GTC", "V"),
    ("GTT", "V"),
    ("TTT", "F"),
    ("TTC", "F"),
    ("TGG", "W"),
    ("TAC", "Y"),
    ("TAT", "Y"),
    ("AAC", "N"),
    ("AAT", "N"),
    ("TGC", "C"),
    ("TGT", "C"),
    ("CAG", "Q"),
    ("CAA", "Q"),
    ("TCT", "S"),
    ("TCC", "S"),
    ("TCA", "S"),
    ("TCG", "S"),
    ("AGC", "S"),
    ("AGT", "S"),
    ("ACT", "T"),
    ("ACC", "T"),
    ("ACA", "T"),
    ("ACG", "T"),
    ("GAT", "D"),
    ("GAC", "D"),
    ("GAA", "E"),
    ("GAG", "E"),
    ("CGT", "R"),
    ("CGC", "R"),
    ("CGA", "R"),
    ("CGG", "R"),
    ("AGG", "R"),
    ("AGA", "R"),
    ("AAA", "K"),
    ("AAG", "K"),
    ("CAT", "H"),
    ("CAC", "H"),
    ("GGT", "G"),
    ("GGC", "G"),
    ("GGA", "G"),
    ("GGG", "G"),
    ("CCT", "P"),
    ("CCC", "P"),
    ("CCA", "P"),
    ("CCG", "P"),
]


CODON_TABLE = {
    "TTT": "F",
    "TCT": "S",
    "TAT": "Y",
    "TGT": "C",
    "TTC": "F",
    "TCC": "S",
    "TAC": "Y",
    "TGC": "C",
    "TTA": "L",
    "TCA": "S",
    "TAA": "_",
    "TGA": "_",
    "TTG": "L",
    "TCG": "S",
    "TAG": "_",
    "TGG": "W",
    "CTT": "L",
    "CCT": "P",
    "CAT": "H",
    "CGT": "R",
    "CTC": "L",
    "CCC": "P",
    "CAC": "H",
    "CGC": "R",
    "CTA": "L",
    "CCA": "P",
    "CAA": "Q",
    "CGA": "R",
    "CTG": "L",
    "CCG": "P",
    "CAG": "Q",
    "CGG": "R",
    "ATT": "I",
    "ACT": "T",
    "AAT": "N",
    "AGT": "S",
    "ATC": "I",
    "ACC": "T",
    "AAC": "N",
    "AGC": "S",
    "ATA": "I",
    "ACA": "T",
    "AAA": "K",
    "AGA": "R",
    "ATG": "M",
    "ACG": "T",
    "AAG": "K",
    "AGG": "R",
    "GTT": "V",
    "GCT": "A",
    "GAT": "D",
    "GGT": "G",
    "GTC": "V",
    "GCC": "A",
    "GAC": "D",
    "GGC": "G",
    "GTA": "V",
    "GCA": "A",
    "GAA": "E",
    "GGA": "G",
    "GTG": "V",
    "GCG": "A",
    "GAG": "E",
    "GGG": "G",
}

CODON_LIST_NNN = CODON_TABLE.values()
CODON_LIST_NNS = [
    "AAC",
    "AAG",
    "ATC",
    "ATG",
    "ACC",
    "ACG",
    "AGC",
    "AGG",
    "TAC",
    "TAG",
    "TTC",
    "TTG",
    "TCC",
    "TCG",
    "TGC",
    "TGG",
    "CAC",
    "CAG",
    "CTC",
    "CTG",
    "CCC",
    "CCG",
    "CGC",
    "CGG",
    "GAC",
    "GAG",
    "GTC",
    "GTG",
    "GCC",
    "GCG",
    "GGC",
    "GGG",
]


# coding DNA sequence of toxin or antitoxin
WT_PARE_DNA = "ATGGCCGTTAGGCTTGTCTGGTCACCCACAGCCAAGGCGGATCTCATCGACATCTACGTGATGATAGGCAGCGAAAACATACGGGCGGCCGACCGCTACTACGATCAGTTGGAAGCGAGGGCCTTACAGCTGGCAGACCAGCCGCGCATGGGTGTCAGGCGACCGGATATCAGGCCTTCCGCACGAATGCTGGTGGAGGCGCCGTTTGTACTGCTCTACGAAACGGTACCGGACACCGACGATGGCCCTGTTGAGTGGGTTGAAATCGTTCGTGTGGTGGATGGGCGGCGGGATCTGAACCGTCTGTTCTAG"
WT_PARD_DNA = "ATGGCAAACGTGGAAAAAATGAGCGTGGCCGTGACCCCGCAACAGGCCGCAGTGATGCGAGAGGCGGTGGAAGCGGGCGAGTATGCCACCGCAAGCGAGATTGTGCGCGAAGCGGTGCGGGATTGGCTGGCCAAGCGCGAACTGCGGCATGACGATATCCGCCGGCTGAGGCAGCTCTGGGATGAAGGCAAAGCAAGCGGGAGACCGGAGCCCGTGGATTTCGACGCGTTGCGAAAGGAAGCTCGGCAAAAGCTGACGGAAGTCCCGCCGAATGGCCGTTAG"

#amino acid sequence of toxin and antitoxin
WT_PARE_DNA_AA = "MAVRLVWSPTAKADLIDIYVMIGSENIRAADRYYDQLEARALQLADQPRMGVRRPDIRPSARMLVEAPFVLLYETVPDTDDGPVEWVEIVRVVDGRRDLNRLF*"
WT_PARD_DNA_AA = "MANVEKMSVAVTPQQAAVMREAVEAGEYATASEIVREAVRDWLAKRELRHDDIRRLRQLWDEGKASGRPEPVDFDALRKEARQKLTEVPPNGR*"

# 5' and 3' flanking sequences containing the bicistron construct
PARE_BCD_90BP_5 = "GGGCCCAAGTTCACTTAAAAAGGAGATCAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGCGACGGTCCCGCCTTGAGCTC"
PARE_BCD_90BP_3 = "AAGCTTGGCTGTTTTGGCGGATGAGAGAAGATTTTCAGCCTGATACAGATTAAATCAGAACGCAGAAGCGGTCTGATAAAACAGAATTTG"
PARD_BCD_90BP_5 = "GGGCCCAAGTTCACTTAAAAAGGAGATCAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGCGACGGTCTCTGAATGAGCTC"
PARD_BCD_90BP_3 = "AAGCTTGGCTGTTTTGGCGGATGAGAGAAGAAATTCGTCGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAGAAGGCCATCCTGACGG"

# barcoded toxin.
WT_PARE_DNA_DNA_CCC = WT_PARE_DNA + "CCC"
WT_PARE_DNA_DNA_TGA = WT_PARE_DNA + "TGA"

PARE_ORF_LENGTH = len(WT_PARE_DNA)
PARD_ORF_LENGTH = len(WT_PARD_DNA)
CODON_START = 90
WT_PARE_DNA_WITH_FLANKS = PARE_BCD_90BP_5 + WT_PARE_DNA + PARE_BCD_90BP_3
WT_PARD_DNA_WITH_FLANKS = PARD_BCD_90BP_5 + WT_PARD_DNA + PARD_BCD_90BP_3
