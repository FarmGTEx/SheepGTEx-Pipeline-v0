#!/bin/bash
num=$1

for ase in `cat $num`
do
    filename=`basename $ase` ; echo $filename ; 
    sed 's/_Aligned.sortedByCoord.out.filtered.rmdup.sorted//g' $ase > gene_ae/$filename
done
