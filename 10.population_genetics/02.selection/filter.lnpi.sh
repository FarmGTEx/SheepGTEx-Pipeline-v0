#!/bin/bash
pop1=$1
pop2=$2
win=$3
step=$4

# filter and extract top 1% region
lines=`cat ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter | wc -l`
top_lines=$(($lines/100+1))
csvtk sort -t -k 7:Nr ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter | head -$top_lines > ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01
sed '1d' ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01 | sort -V | bedtools merge -i - -d 1 -c 4,5,6,7 -o sum,max,max,max | sort -k7 -r | sed "1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI_${pop1}\tPI_${pop2}\tln(PI_${pop1}/PI_${pop2})" > ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged

# variant annotation
awk '{print $1"\t"$2"\t"$3"\t0\t0"}' ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged | sed '1d' > ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged.region
/storage/public/home/2021050427/software/annovar/annotate_variation.pl --outfile ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged.region --buildver sheep ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged.region /storage/public/home/2020060185/genome/sheep/reference/annovar
awk 'BEGIN{print "chr\tstart\tend\ttype\tann"}{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged.region.variant_function > ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged.region.variant_function.addheader
csvtk join -t -f "1,2,3;1,2,3" ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged.region.variant_function.addheader > ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter.t0.01.merged.annotation

