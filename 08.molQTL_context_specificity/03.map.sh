#!/bin/bash
group0=All
group1=Europe
group2=Central_and_East_Asia

# linear regression model (tensorQTL)
## 1. permutation
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	mkdir -p ${tis}/results/tensorqtl/permutation
	for group in ${group0} ${group1} ${group2}
	do
	jsub -q gpu -n 1 -gpgpu "1 mig=1" -R "span[hosts=1]" -J v4.ancqtl_permutation_${tis}_${group} \
		-e ${tis}/log/03.permutation.${tis}.${group}.%J.log -o ${tis}/log/03.permutation.${tis}.${group}.%J.log \
		"bash tensorqtl_permutation.sh ${tis} ${group}"
	done
done

## eGene numbers
echo -e 'Tissue\tGroup\tNumber of eGenes\tNumber of tested genes\tProportion of eGenes' > egenes.txt
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
        do
	zcat ${tis}/results/tensorqtl/permutation/${tis}.${group}.cis_qtl_fdr0.05.txt.gz | sed '1d' | awk -v tis=$tis -v group=$group '{if ($NF=="TRUE"){t+=1}else{f+=1}}END{print tis"\t"group"\t"t"\t"t+f"\t"t/(t+f)}'
	done
done >> egenes.txt

## combine results of permutation
zcat Muscle/results/tensorqtl/permutation/Muscle.${group0}.cis_qtl_fdr0.05.txt.gz | head -1 | sed 's/^/tissue\tgroup\t/g' > tensorqtl_permutation.txt
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	for group in ${group0} ${group1} ${group2}
	do
	zcat ${tis}/results/tensorqtl/permutation/${tis}.${group}.cis_qtl_fdr0.05.txt.gz | awk -v tis=$tis -v group=$group '{if (NR!=1) print tis"\t"group"\t"$0}' >> tensorqtl_permutation.txt
	done
done

## 2. nominal
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	mkdir -p ${tis}/results/tensorqtl/nominal
	for group in ${group0} ${group1} ${group2}
	do
	jsub -q gpu -n 1 -gpgpu "1 mig=1" -R "span[hosts=1]" -J v4.ancqtl_nominal_${tis}_${group} \
		-e ${tis}/log/03.nominal.${tis}.${group}.%J.log -o ${tis}/log/03.nominal.${tis}.${group}.%J.log \
		"bash tensorqtl_nominal.sh ${tis} ${group}"
	done
done

### extract significant cis-SNPs per tisue
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
        for group in ${group0} ${group1} ${group2}
        do
        jsub -q normal -n 1 -R "span[hosts=1]" -J v4.ancqtl_sig_${tis}_${group} \
			-e ${tis}/log/03.sig.${tis}.${group}.%J.log -o ${tis}/log/03.sig.${tis}.${group}.%J.log \
			"bash tensorqtl_sig.sh ${tis} ${group} 1"
	done
done

## 3. independent
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	mkdir -p ${tis}/results/tensorqtl/independent
	for group in ${group0} ${group1} ${group2}
	do
	jsub -q gpu -n 1 -gpgpu "1 mig=1" -R "span[hosts=1]" -J v4.ancqtl_independent_${tis}_${group} \
		-e ${tis}/log/03.independent.${tis}.${group}.%J.log -o ${tis}/log/03.independent.${tis}.${group}.%J.log \
		"bash tensorqtl_independent.sh ${tis} ${group}"
	done
done
### combine results
zcat Muscle/results/tensorqtl/independent/Muscle.${group0}.cis_independent_qtl.txt.gz | head -1 | sed 's/^/tissue\tgroup\t/g' > tensorqtl_independent.txt
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
        for group in ${group0} ${group1} ${group2}
        do
        zcat ${tis}/results/tensorqtl/independent/${tis}.${group}.cis_independent_qtl.txt.gz | awk -v tis=$tis -v group=$group '{if (NR!=1) print tis"\t"group"\t"$0}' >> tensorqtl_independent.txt
        done
done
