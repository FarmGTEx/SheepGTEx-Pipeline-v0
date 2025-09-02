#!/bin/bash
strong_file="strong_pairs.MashR_input.txt.gz"
random_file="random_pairs.MashR_input.txt.gz"

mkdir -p output
jsub -q fat -n 4 -R "span[hosts=1]" -J mashr0 -e log/03.mashr0.%J.log -o log/03.mashr0.%J.log \
    "Rscript run_MashR.R ${strong_file} ${random_file} 0 ./output/top_pairs"
jsub -q fat -n 4 -R "span[hosts=1]" -J mashr1 -e log/03.mashr1.%J.log -o log/03.mashr1.%J.log \
    "Rscript run_MashR.R ${strong_file} ${random_file} 1 ./output/top_pairs_across_all"
