#!/bin/bash
sed '1d' ../tissue40.list | while read tis size num
do
	if [ $size -lt 200 ]; then
		/usr/bin/cp ${tis}/covFile/${tis}.elbow_genoPC5.tsv ${tis}/covFile/${tis}.tsv
	else
		/usr/bin/cp ${tis}/covFile/${tis}.elbow_genoPC10.tsv ${tis}/covFile/${tis}.tsv
	fi
done

# linear regression model (tensorQTL)
## 1. permutation
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/results/tensorqtl/permutation
	jsub -q gpu -n 1 -gpgpu "1 mig=1" -R "span[hosts=1]" -J tensorqtl_permutation_${tis} \
		-e ${tis}/log/04.permutation.${tis}.%J.log -o ${tis}/log/04.permutation.${tis}.%J.log \
		"bash tensorqtl_permutation.sh ${tis}"
done

### eGene numbers
echo -e 'Tissue\tNumber of eGenes\tNumber of tested genes\tProportion of eGenes' > egenes.txt
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	zcat ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | sed '1d' | awk -v tis=$tis '{if ($NF=="TRUE"){t+=1}else{f+=1}}END{print tis"\t"t"\t"t+f"\t"t/(t+f)}'
done >> egenes.txt
awk '$NF=="TRUE"{print $2}' tensorqtl_permutation.txt | sort -u > egene.list
cut -f2 tensorqtl_permutation.txt | sed '1d' | sort -u > allgene.list
sort allgene.list egene.list | uniq -u > nonegene.list

### combine results of permutation
zcat Heart/results/tensorqtl/permutation/Heart.cis_qtl_fdr0.05.txt.gz | head -1 | sed 's/^/tissue\t/g' > tensorqtl_permutation.txt
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	zcat ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | awk -v tis=$tis '{if (NR!=1) print tis"\t"$0}' >> tensorqtl_permutation.txt
done

## 2. nominal
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/results/tensorqtl/nominal
	jsub -q normal -n 2 -R "span[hosts=1]" -J tensorqtl_nominal_${tis} \
		-e ${tis}/log/04.nominal.${tis}.%J.log -o ${tis}/log/04.nominal.${tis}.%J.log \
		"bash tensorqtl_nominal.sh ${tis}"
done
### extract 20k random cis-SNPs per tissue
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J tensorqtl_stat_${tis} \
			-e ${tis}/log/04.stat.${tis}.%J.log -o ${tis}/log/04.stat.${tis}.%J.log \
			"bash tensorqtl_stat.sh ${tis}"
done
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        awk -v tis=$tis '{print tis"\t"$0}' ${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.shuf20k.txt
done | sed '1itissue\tphenotype_id\tvariant_id\tstart_distance\tpval_nominal' > tensorqtl.cis_qtl_pairs.shuf1M.txt

## 3. independent
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
	mkdir -p ${tis}/results/tensorqtl/independent
	jsub -q gpu -n 1 -gpgpu "1 mig=1" -R "span[hosts=1]" -J tensorqtl_independent_${tis} \
		-e ${tis}/log/04.independent.${tis}.%J.log -o ${tis}/log/04.independent.${tis}.%J.log \
		"bash tensorqtl_independent.sh ${tis}"
done
### combine results
zcat Heart/results/tensorqtl/independent/Heart.cis_independent_qtl.txt.gz | head -1 | awk '{print "tissue\t"$0}' | gzip -c > tensorqtl_independent_qtl.txt.gz
for tis in `cut -f1 ../tissue40.list | sed '1d'` ; do zcat ${tis}/results/tensorqtl/independent/${tis}.cis_independent_qtl.txt.gz | awk -v tis=$tis 'NR!=1{print tis"\t"$0}' | gzip -c >> tensorqtl_independent_qtl.txt.gz ; done

### extract significant cis-SNPs per tisue
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J tensorqtl_sig_${tis} \
			-e ${tis}/log/04.sig.${tis}.%J.log -o ${tis}/log/04.sig.${tis}.%J.log \
			"bash tensorqtl_sig.sh ${tis} 1"
done
### combine results
zcat Heart/results/tensorqtl/nominal/Heart.cis_qtl_pairs.sig.txt.gz | head -1 | awk '{print "tissue\t"$0}' | gzip -c > tensorqtl_cis_qtl_pairs.sig.txt.gz
for tis in `cut -f1 ../tissue40.list | sed '1d'` ; do zcat ${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | awk -v tis=$tis 'NR!=1{print tis"\t"$0}' | gzip -c >> tensorqtl_cis_qtl_pairs.sig.txt.gz ; done

## calculate the Inflation Factor (λ)
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        mkdir -p ${tis}/results/tensorqtl/trans
        jsub -q gpu -n 2 -gpgpu "1 mig=1" -R "span[hosts=1]" -J tensorqtl_lambda_${tis} \
			-e ${tis}/log/04.lambda.${tis}.%J.log -o ${tis}/log/04.lambda.${tis}.%J.log \
			"bash tensorqtl_lambda.sh ${tis}"
done
# combine results
echo -e "Tissue\tphenotype_id\tlambdap\tlambdaz" > tensorqtl_lambda.txt
for tis in `cut -f1 ../tissue40.list | sed '1d'` ; do awk -v tis=$tis '{print tis"\t"$1"\t"$2"\t"$3}' ${tis}/results/tensorqtl/trans/${tis}.lambda.txt ; done >> tensorqtl_lambda.txt
