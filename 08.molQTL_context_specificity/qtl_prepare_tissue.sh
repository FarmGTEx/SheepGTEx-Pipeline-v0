#!/bin/bash
rm -f tissue.filtered1.list
for tis in `cut -f1 ../../../tissue40.list | sed '1d'`
do
	grep -w $tis ../../../ind_tis_sample_melt40.tsv | csvtk join -t -H -f '1;1' - ../group.list | awk '$5=="Europe"||$5=="Central_and_East_Asia"{print $1"_"$2"\t"$3"\t"$4"\t"$5}' > ${tis}.grouplist
	code=`cut -f4 ${tis}.grouplist | sort | uniq -c | awk '$1<40' | wc -l`
	if [ $code -eq 0 ] ; then
		echo $tis >> tissue.filtered1.list
       		mkdir -p ${tis}/phenotypes ${tis}/genotypes ${tis}/covFile ${tis}/results ${tis}/log
		mv ${tis}.grouplist ${tis}/${tis}.grouplist
	else
		echo Insufficient sample size in $tis
		rm -f ${tis}.grouplist
	fi
done
