# script allows for classifying a sequencing read into mutations and other categories, see below


from sys import argv
import re
from mutTools import hamming, fasta_iter
from constants import (
    WT_PARE_DNA,
    WT_PARD_DNA,
    PARE_BCD_90BP_5,
    PARE_BCD_90BP_3,
    PARD_BCD_90BP_5,
    PARD_BCD_90BP_3,
    WT_PARE_DNA_DNA_CCC,
    WT_PARE_DNA_DNA_TGA,
)


def mapAndClassify300Read(readName, readSeq, template):
    '''
    classifies a paired end read sequence
    :param readName
    :param readSeq
    :param template: 'pard' or 'pare'
    :return: a read being classified as either of

            'noStartFound':     Could not find the start of the coding sequence defined by the restriction site GAGCTC
                                (ATG should follow right after that)
            'readShort':        the read is too short, shorter than the length of the coding DNA plus flanking DNA
            'wtCCC' or 'wtTGA': some experiments have barcoded wild-type toxin spiked in that are barcoded with an
                                additional 'CCC' or 'TGA' seuqence right after the coding sequence
            'flankMismatch':    whether read sequence flanking the coding sequence has any mutations (indicative of
                                errors in restriction cloning)
            'wt_dna':           wild-type coding sequence
            'hamming>21':       too many mutations
            'n_codons':         n standing for the number of codons mutated, and the mutation string indicative of
                                which nt base is mutated

    '''

    # get right templates
    if template == "pare":
        wt_seq = WT_PARE_DNA
        wt_len = len(wt_seq)
        flank_10bp_5 = PARE_BCD_90BP_5[-10:]
        flank_6bp_3 = PARE_BCD_90BP_3[:6]

    elif template == "pard":
        wt_seq = WT_PARD_DNA
        wt_len = len(wt_seq)
        flank_10bp_5 = PARD_BCD_90BP_5[-10:]
        flank_6bp_3 = PARD_BCD_90BP_3[:6]

    # check where the ATG starts, should be 31, but could also be 30,
    try:
        startPos = re.search("GAGCTC", readSeq[10:]).start() + 16
    except AttributeError:
        toWrite = [readName, "noStartFound"]
        return toWrite

    # find the read length, min read length for fl is 312
    if template == "pare":
        min_len_read = len(wt_seq) + startPos + 13
    elif template == "pard":
        min_len_read = len(wt_seq) + startPos + 10
    readL = len(readSeq)
    if readL < min_len_read:
        toWrite = [readName, "readShort"]
        return toWrite

    # check CCC or TGA wt barcode, only applicable to toxin seq
    if template == "pare":
        if len(readSeq[startPos: startPos + 315]) != len(WT_PARE_DNA_DNA_CCC):
            print(
                readName,
                readSeq,
                len(readSeq[startPos: startPos + 315]),
                len(WT_PARE_DNA_DNA_CCC),
            )

        if hamming(readSeq[startPos: startPos + 315], WT_PARE_DNA_DNA_CCC) == 0:
            toWrite = [readName, "wtCCC"]
            return toWrite

        if hamming(readSeq[startPos: startPos + 315], WT_PARE_DNA_DNA_TGA) == 0:
            toWrite = [readName, "wtTGA"]
            return toWrite
    # check that flanking sequences match up perfectly.
    if (
            hamming(readSeq[startPos - 10: startPos + 3], flank_10bp_5 + "ATG") != 0
            or
            hamming(readSeq[startPos + wt_len - 3: startPos + wt_len + 6], wt_seq[-3:] + flank_6bp_3, ) != 0
    ):
        toWrite = [readName, "flankMismatch"]
        return toWrite

    # find hamming distance of the block
    hammingDist = hamming(readSeq[startPos: startPos + wt_len], wt_seq)

    # find wild type sequence
    if hammingDist == 0:
        toWrite = [readName, "wt_dna_" + template]
        return toWrite

    # random sequence
    elif hammingDist > 21:
        # too many muts
        toWrite = [readName, "hamming>21"]
        return toWrite
    else:
        # create muts
        muts = {}
        c = 0
        for (wtC, mutC) in zip(wt_seq, readSeq[startPos: startPos + wt_len]):
            if wtC != mutC:
                muts[int(c)] = mutC
                # muts.append(str(+c)+':'+mutC)
            c += 1

        numCodonMuts = len(set([int(i) / 3 for i in muts.keys()]))
        fullMutStr = ",".join(str(n) + ":" + mutC for n, mutC in muts.items())
        toWrite = [
            readName,
            str(numCodonMuts) + "_codons",
            "f",
            str(-startPos),
            str(-startPos + readL),
            fullMutStr,
        ]
        return toWrite


if __name__ == "__main__":
    # testing things.
    template, fin, pout = argv[1:4]

    fName = fin.split("/")[-1][:-6]
    fpout = pout + fName + ".csv"

    readLs = []
    with open(fpout, "w") as fout:
        lCount = 0
        for readTuple in fasta_iter(fin):
            # 1 get read name
            readName = readTuple[0][15:27]
            readSeq = readTuple[1]

            readL = len(readSeq)
            readLs.append(readL)
            toWrite = mapAndClassify300Read(readName, readSeq, template)
            fout.write("\t".join(toWrite) + "\n")

    fout.close()
