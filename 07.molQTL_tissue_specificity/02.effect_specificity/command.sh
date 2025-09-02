#!/bin/bash
# 01. prepare tissue combinations
sed '1d' ../tissue40.list | awk '{ for (i = 1; i <= NR; i++) { if (i < NR) { print $1"\t"lines[i] } } lines[NR] = $1 }' > combinations.list
csvtk cut -t -f 'tissue,phenotype_id,variant_id,pval_nominal,slope,slope_se,is_eGene' ../v1.min40_split/tensorqtl_permutation.txt > tensorqtl_permutation.txt
## join tissue combinations
cat combinations.list | while read tis1 tis2
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J 01.join_${tis1}_${tis2} \
		-e combinations/${tis1}_${tis2}/01.join.${tis1}_${tis2}.%J.log -o combinations/${tis1}_${tis2}/01.join.${tis1}_${tis2}.%J.log \
		"bash 01.prepare_joined_permutation.sh $tis1 $tis2"
done

# 02. calculate PPH4 and lead SNP LD between tissue pairs
## lead SNP LD
cat combinations.list | while read tis1 tis2
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J 02.plink_ld_${tis1}_${tis2} \
		-e combinations/${tis1}_${tis2}/02.plink_ld.${tis1}_${tis2}.%J.log -o combinations/${tis1}_${tis2}/02.plink_ld.${tis1}_${tis2}.%J.log \
		"bash 02.plink_ld.sh $tis1 $tis2"
done

## PPH4
cat combinations.list | while read tis1 tis2
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J 02.coloc_pph4_${tis1}_${tis2} \
		-e combinations/${tis1}_${tis2}/02.coloc_pph4.${tis1}_${tis2}.%J.log -o combinations/${tis1}_${tis2}/02.coloc_pph4.${tis1}_${tis2}.%J.log \
		"bash 02.coloc_pph4.sh $tis1 $tis2"
done

# 03. tissue specific effect size
cat combinations.list | while read tis1 tis2
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J 03.het_eqtl_${tis1}_${tis2} \
		-e combinations/${tis1}_${tis2}/03.oppo_eqtl.${tis1}_${tis2}.%J.log -o combinations/${tis1}_${tis2}/03.oppo_eqtl.${tis1}_${tis2}.%J.log \
		"bash 03.tissue_specific_effect_qtl.sh $tis1 $tis2"
done
awk 'NR==1||FNR>1' combinations/*/oppo_effect_eqtl.txt > oppo_effect_eqtl.txt
