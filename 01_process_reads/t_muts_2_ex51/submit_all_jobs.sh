#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 4-23:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=40000M                         # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like


source activate dd_stan_env
source /n/groups/marks/users/david/github/coevolution_mechanism/source_modules.sh

python3 02_usearch_filter.py
python3 03_split_by_at_index.py
python3 04_classify.py
python3 04b_concat_class_files.py
python3 05_count.py
