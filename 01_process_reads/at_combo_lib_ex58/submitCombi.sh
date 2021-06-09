#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 4-23:59                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30G                         # Memory total in MB (for all cores)
#SBATCH -o /n/groups/marks/users/david/ex58/03_called/errs/%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e /n/groups/marks/users/david/ex58/03_called/errs/%j.err                 # File to which STDERR will be written, including job ID

source activate dd_stan_env
module load gcc/6.2.0
python3 process_one_combi.py $1 $2

