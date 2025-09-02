#!/bin/bash
tis=$1
group=$2
threads=$3

# extract significant eQTLs
zcat ${tis}/results/tensorqtl/permutation/${tis}.${group}.cis_qtl_fdr0.05.txt.gz | awk '{if (NR==1 || $NF=="TRUE") print}' > ${tis}/results/tensorqtl/permutation/${tis}.${group}.cis_qtl_fdr0.05.egenes.txt
csvtk -t -j $threads cut -f 'phenotype_id,pval_nominal_threshold' ${tis}/results/tensorqtl/permutation/${tis}.${group}.cis_qtl_fdr0.05.egenes.txt > ${tis}/results/tensorqtl/nominal/${tis}.${group}.pval_nominal_threshold
zcat ${tis}/results/tensorqtl/nominal/${tis}.${group}.cis_qtl_pairs.chr1.txt.gz | csvtk join -t -j $threads -f 'phenotype_id;phenotype_id' - ${tis}/results/tensorqtl/nominal/${tis}.${group}.pval_nominal_threshold | awk '$7<=$10' | gzip -c > ${tis}/results/tensorqtl/nominal/${tis}.${group}.cis_qtl_pairs.sig.txt.gz
for chr in chr{2..26}
do
	echo $chr
	zcat ${tis}/results/tensorqtl/nominal/${tis}.${group}.cis_qtl_pairs.${chr}.txt.gz | csvtk join -t -j $threads -f 'phenotype_id;phenotype_id' - ${tis}/results/tensorqtl/nominal/${tis}.${group}.pval_nominal_threshold | awk '$7<=$10' | sed '1d' | gzip -c >> ${tis}/results/tensorqtl/nominal/${tis}.${group}.cis_qtl_pairs.sig.txt.gz
done
