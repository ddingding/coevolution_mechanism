#pipeline tools
# set of functions for pre-processing of reads, like merging, quality filtering, splitting of read files.

from os import listdir
from os.path import isfile, join
import libClass as lc
import pickle
import os
from collections import defaultdict

from mutTools import fasta_iter_py3
import mapClassPe as mcp

def merge_paired_reads_vsearch_o2(fastq1, fastq2, fout):
    '''
    calls vsearch and merged paired end reads on the o2 cluster.
    requires a path to wherever vsearch is installed
    for vsearch intallation see: https://github.com/torognes/vsearch
    expects 2 paired end read files.
    '''
    vsearchCmd = (
            "/n/groups/marks/users/david/apps/vsearch/bin/vsearch "
            + "--fastq_mergepairs %s  --reverse %s --fastqout %s_merged.fastq"
            % (fastq1, fastq2, fout)
    )
    print(vsearchCmd)
    os.system(vsearchCmd)
    return vsearchCmd


def filter_fastq_quality(fastq_in, fasta_out):
    """
    Take a path to fastq_in, and specify a path to fasta_out,
    and filter merged fastq files for quality on the o2 cluster
    requires a path to wherever vsearch is installed
    for vsearch intallation see: https://github.com/torognes/vsearch
    """
    print(fastq_in, fasta_out)

    # new style formatting
    vsearchCmd = "/n/groups/marks/users/david/apps/vsearch/bin/vsearch --fastq_filter {0} --fastq_truncqual 20 " \
                 "--fastq_maxns 3 --fastq_maxee 0.5 --fastq_ascii 33 --fastaout {1}.fasta".format(
        fastq_in, fasta_out
    )

    print(vsearchCmd)
    os.system(vsearchCmd)
    return vsearchCmd


def split_by_indices(fa_f, index_to_at, dout, f_name, verbose=True):
    """
    Take a fasta file, and split it by the DNA sequence after the restriction site.

    -------INPUT------------
    fa_f            path to fasta_file to split
    index_to_at     dictionary of str index to look for in sequence
                    after the restriciton site AAGCTT, to str mutkey
    dout            str, directory
    f_name          prefix to add for the filenames of the split files,
                    eg. primer name
    """
    if verbose:
        print("split_by_index for", fa_f, index_to_at, dout, f_name)

    # make a file to write for each index
    ind_to_file = {}
    for i, at_mut in index_to_at.items():
        ind_to_file[i] = open(dout + f_name + "_" + at_mut + ".fasta", "w")
        print("created file ", dout + f_name + "_" + at_mut + ".fasta")
    # to write all the index mismatched sequences to
    ind_to_file["unexpected"] = open(dout + f_name + "_unexpected.fasta", "w")

    print(fa_f, " looking for indices: ", ind_to_file.keys())
    # iter through fasta file and write into appropriate file

    unexpected_indices = defaultdict(int)
    not_found = 0
    total_seqs = 0

    for fa_rec in fasta_iter_py3(fa_f):
        n = fa_rec.header
        s = fa_rec.sequence

        # looking for where the AAGCTT starts in the sequence
        re_site_pos = s.find("AAGCTT")
        # if re site not found in string
        if re_site_pos == -1:
            at_ind = str(None)
            not_found += 1

        else:
            at_ind = s[re_site_pos + 6: re_site_pos + 10]

        # check whether the index is expected
        if at_ind in ind_to_file.keys():
            ind_to_file[at_ind].write(">" + n + "_" + at_ind + "\n" + s + "\n")
        else:
            unexpected_indices[at_ind] += 1
            ind_to_file["unexpected"].write(">" + n + "_" + at_ind + "\n" + s + "\n")

        total_seqs += 1

    # write the summary files
    with open(dout + f_name + "_split_stats.csv", "w") as statOut:
        print(unexpected_indices)
        print(sum(unexpected_indices.values()))
        statOut.write(
            "unexpected indices not split: "
            + str(sum(unexpected_indices.values()))
            + "\n"
        )
        statOut.write(
            "# seq split: " + str(total_seqs - sum(unexpected_indices.values())) + "\n"
        )

        # sorted
        sorted_unexpected_indices = sorted(
            unexpected_indices.items(), key=lambda kv: kv[1], reverse=True
        )
        statOut.write("unexpected indices are: " + str(sorted_unexpected_indices))
    # close files
    for f in ind_to_file.values():
        f.close()


def demultiplex_fastas(list_fastas, fasta_dout, df_config, config_dics, exp_num=1):
    '''
    Samples of toxin single mutants in the background of different antitoxin mutants were pooled in the same flask.
    The antitoxin mutant background of a particular toxin mutant is encoded on the barcode that is 3' to the toxin
    stop codon.
    This function allows for splitting a list of fasta files based on that barcode.
    For a given fasta file, it references from the
    filename (BioMicroCenter (BMC) name) --> my primer number --> which cell mix --> at_indices --> AT background
    and then calls the split_by_indices()
    :param list_fastas: absolute path fasta files to demultiplex, such as
    '/n/groups/marks/users/david/ex47/03filtered_2_3786W/190311Lau_D19-2187_merged.fasta'
    :param fasta_dout: directory to write files out
    :param exp_num:
    :param df_config:
    :param config_dics:
    :return:
    '''

    # primer nums that are toxins
    toxin_primer_nums = list(df_config.loc[df_config.gene == "t"].primer)

    for fa_f in list_fastas:
        # go from filename (BMC name) --> my primer --> which cell mix --> at_indices
        print("demultiplexing", fa_f)
        f_name = fa_f.split("/")[-1]
        # for x47 file naming
        if exp_num  == 1:
            bmc_index = f_name.split('_')[1]
        # for x51 file naming
        else:
            bmc_index = f_name.split("_")[-1][:-2]


        primer_str = map_bmc_to_primer(config_dics, bmc_index)
        # check it is a toxin and needs to be demultiplexed
        if float(primer_str) in toxin_primer_nums:
            # fetch cell_mix number by my primer number
            cell_mix = map_primer_to_cell_mix(df_config, primer_str)
            print(primer_str, cell_mix)
            if exp_num == 1:
                at_index_to_mutkey = map_cell_mix_to_at_indices_1(config_dics, cell_mix)
            else:
                at_index_to_mutkey = map_cell_mix_to_at_indices_2(config_dics, cell_mix)

            print(at_index_to_mutkey)
            if at_index_to_mutkey:
                f_name = str(primer_str)
                # for x51
                fa_f = fa_f + "_.extendedFrags.fasta"
                split_by_indices(fa_f, at_index_to_mutkey, fasta_dout, f_name)


def map_bmc_to_primer(config_dics, bmc_index):
    # expects bmc_index like D19-2181
    # returns my own primer name like 187
    bmc_index_w_lane_1_num = bmc_index + "-1"

    # only first 3 characters are primer, the rest is the index sequence
    primer_num = config_dics["BMC_TO_PRIMER"][bmc_index_w_lane_1_num][:3]
    return primer_num

def map_primer_to_cell_mix(df_config, primer_str):
    # expect a primer number like 178
    # and returns whether this sample has cell mix 1, or cell mix 2
    return int(df_config.loc[df_config.primer == int(primer_str), "cell_mix"].iloc[0])

def map_cell_mix_to_at_indices_1(config_dics, cell_mix):
    """this is for ex47
    # cell_mix is either 1, 2, or 249 (the toxin sample that is not a pool,
    # but just the top10 miniprep)
    """

    if cell_mix_str not in map(str,[1,2,249]):
        print('supplied cell_mix', str(cell_mix), 'doesnt have a mix of at indices.')
        return None

    if cell_mix_str == '1':
        #if all mutkeys present
        at_index_to_mutkey = config_dics['AT_INDEX_TO_MUTKEY']
        return at_index_to_mutkey

    elif cell_mix_str == '2':
        #if just subset of mix 2 is present.
        indices_to_use = config_dics['CELL_MIX_TO_AT_INDICES'][cell_mix]
        at_index_to_mutkey = dict([(k,v) for k,v in
                config_dics['AT_INDEX_TO_MUTKEY'].items() if k in indices_to_use])
        return at_index_to_mutkey

    elif cell_mix_str == '249':
        return None
def map_cell_mix_to_at_indices_2(config_dics, cell_mix):
    # this is for ex51
    # cell mix is either of 1,2,3,4, and -1 for -AT samples that were not pooled in the same flask
    # returns at_indices to split by for that sample
    cell_mix_str = str(cell_mix)

    # this is for ex51
    if cell_mix_str not in map(str, [-1, 1, 2, 3, 4]):
        return None
    elif cell_mix_str == "-1":
        # to catch mcs sample that was not multiplexed
        return None
    else:
        at_index_to_mutkey = config_dics["AT_INDEX_TO_MUTKEY_" + cell_mix_str]
        return at_index_to_mutkey

#### for 04_classify.py
def classify_fasta(fin, fout, template):
    """
    expects:
    fin         fpath to fasta file
    fout        fpath to classify file out
    template    either 'pare' or 'pard'
    """

    print("classifying..." + fin)
    c = 0

    with open(fout, "w") as fout1:
        for fa_rec in fasta_iter_py3(fin):
            n = fa_rec.header
            s = fa_rec.sequence
            writeL = mcp.mapAndClassify300Read(n, s, template)
            fout1.write("\t".join(writeL) + "\n")
            c += 1

            if c % 10000 == 1:
                print(c)


def map_primer_to_template(primer_str, df_config):
    # get primer, 178, and return 'pare'
    t_or_at = df_config.loc[df_config.primer == int(primer_str), "gene"].iloc[0]

    if t_or_at == "t":
        return "pare"
    elif t_or_at == "at":
        return "pard"
    else:
        print("couldnt parse: ", t_or_at)


def class_samples(fasta_dir_in, dout, df_config):
    # take a fasta_dir_in and classifies all .fasta files in that directory into dout.
    fastas_all = [
        f
        for f in listdir(fasta_dir_in)
        if isfile(join(fasta_dir_in, f))
           and f.endswith(".fasta")
           and not f.endswith("unexpected.fasta")
    ]

    done_fs = [
        f.rstrip("_class.tsv")
        for f in listdir(dout)
        if isfile(join(dout, f)) and f.endswith("_class.tsv")
    ]

    fastas = [f for f in fastas_all if f.rstrip(".fasta") not in done_fs]

    print(fastas_all, done_fs, fastas)
    print(len(fastas_all), len(done_fs), len(fastas))
    print("fastas to be done", len(fastas), fastas)
    """
    fastas1 = ['188_K63D.fasta', '184_A16K.fasta',..., '193_K63D.fasta']
    """
    for f in fastas:
        template = map_primer_to_template(f[:3], df_config)
        fin = fasta_dir_in + f
        fout = dout + f[:-6] + "_class.tsv"
        classify_fasta(fin, fout, template)


def concat_fs_in_dirs(dir1, dir2, dout):
    f1s = [f for f in listdir(dir1) if isfile(join(dir1, f))]
    f2s = [f for f in listdir(dir2) if isfile(join(dir2, f))]
    intersect_files = set(f1s).intersection(f2s)
    print("concatenating files:", intersect_files)

    c = 0
    for f in intersect_files:
        print(c, "/", len(intersect_files), "files.")
        concat_files([dir1 + f, dir2 + f], dout + f)
        c += 1
    return


def concat_files(list_files, f_path_out):
    # expect list of files, and merges the files into the f_path_out
    if list_files:
        print("concatenating files:", list_files)
        with open(f_path_out, "w") as fout:
            for f in list_files:
                with open(f, "r") as fin:
                    for l in fin:
                        fout.write(l)

    else:
        print("Error: merge_2_files() wasnt supplied with file list.")


# for 05 make sample objects


def map_primer_to_time(primer_str, df_config):
    t = df_config.loc[df_config.primer == int(primer_str), "t"].iloc[0]
    return t


def fetch_sample_obj_args(classIn, df_config):
    # get a filepath to the classified file, and construct based on the df_config
    # all the parameters needed for making a sample object.
    primer_str = classIn.split("/")[-1].split("_")[0]

    timepoint = map_primer_to_time(primer_str, df_config)
    # make a sample name that has the class file name and the timepoint
    sample_n = classIn.split("/")[-1].rstrip("_class.tsv") + "_t" + str(timepoint)
    template = map_primer_to_template(primer_str, df_config)
    od = 0
    return timepoint, sample_n, template, od


def make_sample_obj(classIn, pickle_out, df_config):
    # can also pass classIn as a list of filepaths to SampleObject to merge
    # into one sample object
    # classIn should look like /n/groups/marks/users/david/ex47/05_class/187_mcsAT_class.tsv

    timepoint, sample_n, template, od = fetch_sample_obj_args(classIn, df_config)

    currObj = lc.SampleObject(
        [classIn],
        sample_n=sample_n,
        timepoint=int(timepoint),
        template=template,
        od=od,
        use_nns=False,
    )
    pickle.dump(currObj, open(pickle_out, "wb"))

    return currObj


def make_all_sample_obj(din, pickle_out, df_config):
    # make sample objects from class.tsv file
    # and save as pickles
    dic_sample_obj = {}
    class_fs = [
        f for f in listdir(din) if isfile(join(din, f)) and f.endswith("class.tsv")
    ]
    # one was missing
    class_fs = ["572_wtAT_class.tsv", "574_F73A_class.tsv"]
    for f in class_fs:
        s_name = f.rstrip("_class.tsv")  # should be 178_wtAT
        fout = pickle_out + s_name + "_sample_obj.p"
        if not isfile(fout):
            print("creating sample", s_name)
            s_obj = make_sample_obj(din + f, fout, df_config)
            dic_sample_obj[s_name] = s_obj
        else:
            print("skipped file", f)
    return dic_sample_obj


# 2 make timecourse obects.
# save as pickles
def make_tc_obj(s_obj1, s_obj2, template, pickle_out_dir):
    tc_name = "_".join(s_obj1.sample_n.split("_")[:-1])

    if not isfile(pickle_out_dir + tc_name + "_tc.p"):
        tc = lc.Timecourse(template, [s_obj1, s_obj2])

        tc.calculateFitness(first_timepoint=s_obj1.timepoint, fit_wrt_sample_stop=True)
        tc.samples[-1].calc_average_aa_fit()

        # for the tc name, just take the sample name and
        # remove the t0 at end
        pickle.dump(tc, open(pickle_out_dir + tc_name + "_tc.p", "wb"))
        print("just pickled dumped", tc_name, "yay!!")

        return tc
    else:
        print("skipping: " + tc_name)


def map_primer_to_gene_sample(primer_str, df_config):
    g = df_config.loc[df_config.primer == int(primer_str), "gene"].iloc[0]
    sample = df_config.loc[df_config.primer == int(primer_str), "sample"].iloc[0]
    return str(g) + "_" + str(sample)


def make_all_possible_tcs_single(sample_list, pickle_out_dir, df_config):
    # expects sample_list := ['547_W59L',...]
    tc_list = []

    # dic_obj is a dictionary from s_name to sample_obj
    # find the objects that are part of the same sample, gene, and mutant
    dic_id_to_s_names = {}
    for s_name in sample_list:
        print(s_name)
        primer = s_name.split("_")[0]
        if len(s_name) > 3:
            mut = s_name.split("_")[1]
        else:
            mut = "m"
            s_name = s_name + "_m"
        gene_sample = map_primer_to_gene_sample(primer, df_config)
        id = gene_sample + "_" + mut

        try:
            dic_id_to_s_names[id].append(s_name)
        except KeyError:
            dic_id_to_s_names[id] = [s_name]

    print("dic_id_to_s_names", dic_id_to_s_names)

    c = 0
    for id, s_names in dic_id_to_s_names.items():

        c += 1
        if isfile(pickle_out_dir + s_names[0] + "_tc.p") or isfile(
                pickle_out_dir + s_names[1] + "_tc.p"
        ):
            print("skipping s_names", c)
            continue

        # construct dic_samples: s_name to sample_obj

        dic_samples = {}
        for s_n in s_names:
            print(s_n)

            dic_samples[s_n] = pickle.load(
                open(pickle_out_dir + s_n + "_sample_obj.p", "rb")
            )

        # make the tc
        if len(s_names) == 2:
            s1_name = min(s_names, key=lambda x: dic_samples[x].timepoint)
            s2_name = max(s_names, key=lambda x: dic_samples[x].timepoint)

            s1_obj = dic_samples[s1_name]
            s2_obj = dic_samples[s2_name]

            template = map_primer_to_template(s_names[0][:3], df_config)

            print("trying to make", s1_name, " ", c, "tc out of", len(sample_list) / 2)
            tc = make_tc_obj(s1_obj, s2_obj, template, pickle_out_dir)
            # tc_list.append(tc)

    return tc_list



def make_all_possible_tcs(dic_samples, pickle_out_dir, df_config):
    #expects dic_samples := 547_W59L : sample_obj
    tc_list = []

    #dic_obj is a dictionary from s_name to sample_obj
    #find the objects that are part of the same sample, gene, and mutant
    dic_id_to_s_names = {}
    for s_name, s_obj in dic_samples.items():
        primer, mut = s_name.split('_')
        gene_sample = map_primer_to_gene_sample(primer, df_config)
        id = gene_sample + '_'+ mut

        try:
            dic_id_to_s_names[id].append(s_name)
        except KeyError:
            dic_id_to_s_names[id] =[s_name]

    print(dic_id_to_s_names)

    for id,s_names in dic_id_to_s_names.items():

        if len(s_names) ==2:
            s1_name = min(s_names, key=lambda x: dic_samples[x].timepoint)
            s2_name = max(s_names, key=lambda x: dic_samples[x].timepoint)
            s1_obj = dic_samples[s1_name]
            s2_obj = dic_samples[s2_name]

            template = map_primer_to_template(s_names[0][:3], df_config)

            tc = make_tc_obj(s1_obj, s2_obj, template, pickle_out_dir)
            tc_list.append(tc)
    return tc_list


def pair_replicate_tc_single(tc_name_list, df_config):
    df_rep = df_config[["primer", "replicate"]].dropna()

    df_rep_dic = dict(zip(df_rep.primer.astype(int), df_rep.replicate.astype(int)))
    print(df_rep_dic)
    # df_rep_dic_2 = dict(zip(df_rep.replicate.astype(str), df_rep.primer.astype(str)))

    # take a list of tc_names
    tc_name_paired_list = []
    # pair them
    for tc_name in tc_name_list:
        tc_pri = int(tc_name[:3])
        # rep_pri =  str(int(df_rep.loc[df_rep['primer'] == int(tc_pri), 'replicate'].iloc[0]))
        # rep_pri =  df_rep.loc[df_rep['primer'] == int(tc_pri), 'replicate'].iloc[0]#['replicate']##.replicate[0]
        if tc_pri in df_rep_dic:
            rep_pri = df_rep_dic[tc_pri]
            print(rep_pri)
            rep_tc_name = str(rep_pri) + tc_name[3:]

            tc_name_paired_list.append((tc_name, rep_tc_name))

        else:
            print("skipping", tc_pri)

    # return list of paired tcs
    return tc_name_paired_list


def pair_replicate_tc(tc_list, df_config):
    # takes a list of timecourse objects
    # returns a list of paired up timecourses
    #
    # take the first primer, which is in the sample_n of s_obj1

    # go through the tcs.

    # df going from primer to replicate, dropped rows with nan values
    # goes from primer number to replicate primer number
    df_rep = df_config[["primer", "replicate"]].dropna()

    def make_s1_name_to_tc(tc_list):
        # make a dictionary with keys {178_wtAT: TimeCourseObject}
        s1_name_to_tc = {}
        for tc in tc_list:
            s_obj1 = tc.samples[0]
            s1_name = s_obj1.sample_n
            s1_name_to_tc[s1_name] = tc

        return s1_name_to_tc

    #
    s1_name_to_tc = make_s1_name_to_tc(tc_list)

    tc_paired_list = []
    for tc in tc_list:
        s1_name = tc.samples[0].sample_n
        s1_pri = s1_name[:3]
        if int(s1_pri) in df_rep.primer.values:
            rep_s1_pri = str(
                int(df_rep.loc[df_rep["primer"] == int(s1_pri), "replicate"].iloc[0])
            )
            rep_s1_name = rep_s1_pri + s1_name[3:]

            tc_paired_list.append((s1_name_to_tc[s1_name], s1_name_to_tc[rep_s1_name]))

    return tc_paired_list

