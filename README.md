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
5) count the mutants, create a table detailing pre- and post selection read counts

Steps:
1) need to set python path for loading modules found in ./src/, so run: source source_modules.sh
2) install flash and vsearch and set the path in ./src/constants.py
3) download files from sra: https://www.ncbi.nlm.nih.gov/sra/PRJNA736482, unzip all the files
3b) for antitoxin singles mutant files only: you will have to deinterleave the .fastq files, which contain paired end reads, into separate files that each contain one of the paired end reads.
for example for sample '180716Lau_D18-6083_phiX_bestmap.fastq', run: 
paste - - - - - - - - < 180716Lau_D18-6083_phiX_bestmap.fastq | tee >(cut -f 1-4 | tr "\t" "\n" > 180716Lau_D18-6083_phiX_bestmap_1_sequence.fastq) |  cut -f 5-8 | tr "\t" "\n" > 180716Lau_D18-6083_phiX_bestmap_2_sequence.fastq
4) set working directory (where all the processed files should land) and download directory in ./src/constants.py
5) run ./setup_dirs.py to create directory structure and move files to the right place.
6) run the scripts in each folder of ./01_process_reads/

The output of these scripts is supplied in ./01_process_reads/output/, and can then be used to perform Bayesian inference of relative growth rates for each amino acid variant.

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
