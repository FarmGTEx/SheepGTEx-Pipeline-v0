#!/bin/bash
sed '1d' tissue100.list | while read tis size num
do
	for group in group1 group2
	do
		if [ $size -lt 400 ]; then
		/usr/bin/cp ${tis}/${group}/covFile/${tis}.elbow_genoPC5.tsv ${tis}/${group}/covFile/${tis}.tsv
		else
		/usr/bin/cp ${tis}/${group}/covFile/${tis}.elbow_genoPC10.tsv ${tis}/${group}/covFile/${tis}.tsv
		fi
	done
done

# linear regression model (tensorQTL)
## 1. permutation
for tis in `cut -f1 tissue100.list | sed '1d'`
do
	for group in group1 group2
	do
	mkdir -p ${tis}/${group}/results/tensorqtl/permutation
	jsub -q normal -n 1 -R "span[hosts=1]" -J tensorqtl_permutation_${tis}_${group} \
		-e ${tis}/${group}/log/03.permutation.${tis}.%J.log -o ${tis}/${group}/log/03.permutation.${tis}.%J.log \
		"bash tensorqtl_permutation.sh ${tis} ${group}"
	done
done

## 2. nominal
for tis in `cut -f1 tissue100.list | sed '1d'`
do
	for group in group1 group2
	do
	mkdir -p ${tis}/${group}/results/tensorqtl/nominal
	jsub -q normal -n 1 -R "span[hosts=1]" -J tensorqtl_nominal_${tis}_${group} \
		-e ${tis}/${group}/log/03.nominal.${tis}.%J.log -o ${tis}/${group}/log/03.nominal.${tis}.%J.log \
		"bash tensorqtl_nominal.sh ${tis} ${group}"
	done
done
