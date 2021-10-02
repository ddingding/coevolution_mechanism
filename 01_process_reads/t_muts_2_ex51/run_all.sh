#!/bin/bash

source activate dd_stan_env
source /n/groups/marks/users/david/github/coevolution_mechanism/source_modules.sh

echo 'running pipeline all'
python3 01_usearch_merge.py
python3 02_usearch_filter.py
python3 03_split_by_at_index.py
python3 04_classify.py
python3 04_classify_no_at.py
python3 04b_concat_class_files.py
python3 04b_concat_class_files_no_at.py
python3 05_count.py
python3 05_count_no_at.py



