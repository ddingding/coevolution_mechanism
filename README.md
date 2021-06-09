# coevolution_mechanism
code for paper

Disclaimer: This code is for reproducibility purposes, and contains multiple legacy features. 

Steps
1: process raw reads to get read counts per codon mutant for each sample
2: run bayesian inference to infer growth rates of amino acid mutants for AT single mutants and T mutants.
3: fit a nonlinear global epistasis model to the data.
4: visualize and further analysis of the data.

For 1): see the analysis_file_pairings.xlsx for the pre-processing code that goes with a particular file.
This code as is will not run on your local machine as it is tied to our particular architecture.

