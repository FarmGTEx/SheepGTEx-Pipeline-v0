#!/bin/bash
pop1=$1
pop2=$2
win=$3
step=$4

#!/bin/bash
mkdir -p ${pop1}_${pop2}
# calculate lnPi ratio
csvtk join -t -f "1,2,3;1,2,3" ${pop1}/chrAuto.${win}_${step}.windowed.pi ${pop2}/chrAuto.${win}_${step}.windowed.pi | sed '1d' | awk '{{if($4>=$6) min=$6; else min=$4} ; printf $1"\t"$2"\t"$3"\t"min"\t"$5"\t"$7"\t%.6f\n",log($5/$7)}' | sed "1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI_${pop1}\tPI_${pop2}\tln(PI_${pop1}/PI_${pop2})" > ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi
awk '$4>=40' ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi > ${pop1}_${pop2}/chrAuto.${win}_${step}.lnPi.filter
