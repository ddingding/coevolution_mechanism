# coevolution_mechanism

This code is for reproducibility purposes to accompany the paper. 

Steps in the processing:
###1. Raw read processing. 
process raw reads to get read counts per codon mutant for each sample. 

The code in each subfolder (./01_process_reads/at_combo_lib/, ./01_process_reads/t_muts_1/, ./01_process_reads/t_muts_2/, ./01_process_reads/at_singles/) goes from raw read files (.fastq) to a Timecourse object.
The processing steps are:
1) merge paired end reads
2) quality filter reads
3) optional: split files by antitoxin barcodes, if toxin single mutants have been run pooled in a flask
4) classify each read into different categories, such as single codon variants
5) convert these classified reads into structured python objects: Sample Object for each sample and timepoint, or Timecourse object, which links Sample objects for a particular flask across timepoints
6) convert these sample or time course objects into tabular form detailing pre- and post selection read counts. ./01_process_reads/06_make_raw_count_df.py 

Step 5 is legacy, and could be skipped by slightly rewriting functions form step 6) to directly count variants on the step 4 output, but it's here for reproducibility purposes.

The output of these scripts is supplied in ./01_process_reads/output/, and can then be used to perform Bayesian inference of relative growth rates for each amino acid variant.

need to set python path for loading modules found in ./src/
install flash and vsearch (for reproducibility reasons)

for example:
download files from sra: https://www.ncbi.nlm.nih.gov/sra/PRJNA736482
unzip all the files
for antitoxin singles mutant files only: you will have to split the .fastq files, which contain paired end reads, into separate paired end read files
set working directory (where all the processed files should land) and download directory in constants.py
run setup_dirs.py to create directory structure and move files to the right place.

See the analysis_file_pairings.xlsx for the pre-processing code that goes with a particular file, or the initial script in each subfolder.


###2. Bayesian inference of growth rates for each amino acid variant.
To go from read counts per codon mutant in each sample, to posterior beliefs of growth rates for amino acid variants.
Briefly, we impart all the hierarchical structure we know about the experiment (such as synonymous codon variants should inform the shared amino acid variant growth rates, and that we should use read count data from both replicate experiments to give us a belief about a particular amino acid variant growth rate) in the model.
Please refer to the paper where we validate this model by comparing our inference on synthetic data, as well as doing extensive posterior predictive checks with our observed data.
Slight tweaks should make this model run on other single mutant datasets, feel free to contact me for that.

###3. Fitting a non-specific, nonlinear model to single and double mutant data. 
This python notebook (can be run on Google Colab) details how we fit a simple independent, nonlinear model the the single and double mutant data as a null model for our expected double mutant growth rates.
The independent part of the model means that we only use independent mutant effect terms for each amino acid variant, and do not try to model epistatic, pairwise terms. You can see in the paper that this simple model does quite well.


For questions, please don't hesitate to contact me via email: firstnamelastname {at} g.harvard.edu
