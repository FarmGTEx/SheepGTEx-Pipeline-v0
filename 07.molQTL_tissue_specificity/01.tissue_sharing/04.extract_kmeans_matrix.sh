#!/bin/bash
jsub -q normal -n 1 -R "span[hosts=1]" -J extract_kmeans_matrix_e \
    -e log/04.extract_kmeans_matrix_e.%J.log -o log/04.extract_kmeans_matrix_e.%J.log \
    "python3 extract_kmeans_matrix.py"
