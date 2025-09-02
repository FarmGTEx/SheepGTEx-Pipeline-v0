#!/bin/bash
pop1=$1
pop2=$2
win=$3
step=$4

# merge overlapping regions
sed '1d' ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01 | sort -V | bedtools merge -i - -d 1 -c 4,5,6 -o sum,max,max | sort -k5 -r | sed '1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST' > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged

# variant annotation
awk '{print $1"\t"$2"\t"$3"\t0\t0"}' ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged | sed '1d' > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged.region
/storage/public/home/2021050427/software/annovar/annotate_variation.pl --outfile ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged.region --buildver sheep ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged.region /storage/public/home/2020060185/genome/sheep/reference/annovar
awk 'BEGIN{print "chr\tstart\tend\ttype\tann"}{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged.region.variant_function > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged.region.variant_function.addheader
csvtk join -t -f "1,2,3;1,2,3" ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged.region.variant_function.addheader > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01.merged.annotation
