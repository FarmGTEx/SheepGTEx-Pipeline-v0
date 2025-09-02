#!/bin/bash
pop1=$1
pop2=$2
win=$3
step=$4

for chr in chr{1..26}
do
    sed '1d' ${pop1}_${pop2}/${chr}.${win}_${step}.windowed.weir.fst | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t%.6f\t%.6f\n",$5,$6}'
done | sed '1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST' > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst
awk '$4>=40' ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter

# filter and extract top 1% region
lines=`cat ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter | wc -l`
top_lines=$(($lines/100+1))
csvtk sort -t -k 5:Nr ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter | \
head -$top_lines > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01
