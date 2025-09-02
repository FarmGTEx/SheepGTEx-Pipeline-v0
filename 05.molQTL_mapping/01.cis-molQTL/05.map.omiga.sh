#!/bin/bash

## 3.2 linear mixed model (omiga)
# cis-QTL mapping
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	jsub -q normal -n 4 -R "span[hosts=1]" -J omiga_cis_${tis} \
	-e ${tis}/log/05.omiga_cis.${tis}.%J.log -o ${tis}/log/05.omiga_cis.${tis}.%J.log \
	"bash omiga_cis.sh ${tis} 4"
done
# combine results
echo -e "Tissue\tphenotype_id\tr_Slope\tr_Z\tr_-log10P\tcis_num\tvariant_id\tpval_g1\teff_g1\tZ\tpval_nominal\tslope_perm\tz_perm\tpval_beta\tis_eGene" > omiga_cis.stat.txt
cat */results/omiga/cis_lmm/omiga.stat.txt >> omiga_cis.stat.txt

### extract significant cis-SNPs per tisue
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J omiga_sig_${tis} \
			-e ${tis}/log/05.sig.${tis}.%J.log -o ${tis}/log/05.sig.${tis}.%J.log "bash omiga_sig.sh ${tis} 1"
done
### combine results
zcat Heart/results/omiga/cis_lmm/Heart.cis_qtl_pairs.sig.txt.gz | head -1 | awk '{print"tissue\t"$0}' | gzip -c > omiga_cis_qtl_pairs.sig.txt.gz
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	zcat ${tis}/results/omiga/cis_lmm/${tis}.cis_qtl_pairs.sig.txt.gz | awk -v tis=$tis 'NR!=1{print tis"\t"$0}' | gzip -c >> omiga_cis_qtl_pairs.sig.txt.gz
done

# calculate lambda
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	jsub -q normal -n 4 -R "span[hosts=1]" -J omiga_lambda_${tis} \
		-e ${tis}/log/05.omiga_lambda.${tis}.%J.log -o ${tis}/log/05.omiga_lambda.${tis}.%J.log \
		"bash omiga_lambda.sh ${tis} 4"
done
echo -e "Tissue\tphenotype_id\tlambdap\tlambdaz" > omiga_lambda.txt
for tis in `cut -f1 ../tissue40.list | sed '1d'` ; do awk -v tis=$tis '{print tis"\t"$1"\t"$2"\t"$2}' ${tis}/results/omiga/lambda/${tis}.GC_lambda.txt | sed '1d' ; done >> omiga_lambda.txt
