#!/bin/bash
k=$1
threads=$2

/storage/public/home/2020060185/software/admixture_linux-1.3.0/admixture -j$threads --cv chrAuto.filtered.keep.pruned.bed $k
