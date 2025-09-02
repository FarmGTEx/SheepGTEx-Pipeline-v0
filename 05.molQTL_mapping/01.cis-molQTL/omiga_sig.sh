#!/bin/bash
tis=$1
threads=$2

# extract significant eQTLs
zcat ${tis}/results/omiga/cis_lmm/${tis}.cis_qtl.txt.gz | awk 'NR==1||$NF<=0.05' | csvtk -t -j $threads cut -f 'pheno_id,pval_g1_threshold' > ${tis}/results/omiga/cis_lmm/${tis}.pval_nominal_threshold
zcat ${tis}/results/omiga/cis_lmm/${tis}.cis_qtl_pairs.chr1.txt.gz | csvtk join -t -j $threads -f 'pheno_id;pheno_id' - ${tis}/results/omiga/cis_lmm/${tis}.pval_nominal_threshold | awk '$7<=$8' | gzip -c > ${tis}/results/omiga/cis_lmm/${tis}.cis_qtl_pairs.sig.txt.gz
for chr in chr{2..26}
do
	echo $chr
	zcat ${tis}/results/omiga/cis_lmm/${tis}.cis_qtl_pairs.${chr}.txt.gz | csvtk join -t -j $threads -f 'pheno_id;pheno_id' - ${tis}/results/omiga/cis_lmm/${tis}.pval_nominal_threshold | awk '$7<=$8' | sed '1d' | gzip -c >> ${tis}/results/omiga/cis_lmm/${tis}.cis_qtl_pairs.sig.txt.gz
done
