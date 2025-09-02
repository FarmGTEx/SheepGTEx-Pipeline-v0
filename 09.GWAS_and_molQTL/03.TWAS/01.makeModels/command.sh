#!/bin/bash
# referred to: https://github.com/hakyimlab/PredictDB-Tutorial
#              https://github.com/FarmGTEx/PigGTEx-Pipeline-v0/blob/master/14_GWAS_and_molQTL/TWAS/TWAS.md

# 1. Preprocess the phenotype,covariate and genotype
## molecular phenotype and covatiate annotation
bash 01.makeAnnotation.sh

# 2. SNP and genotype annotation
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	for qtl in eQTL sQTL eeQTL isoQTL stQTL 3aQTL enQTL
	do
		mkdir -p ${qtl}/${tis}/log
		jsub -q normal -n 1 -R "span[hosts=1]" -J prepare_genotype_${qtl}_${tis} \
			-e ${qtl}/${tis}/log/02.prepare.%J.log -o ${qtl}/${tis}/log/02.prepare.%J.log \
			"bash 02.prepare.sh ${qtl} ${tis}"
	done
done

# 3. Train the models
for tis in `cut -f1 tissue40.list | sed '1d'`
do
        for qtl in eQTL sQTL eeQTL isoQTL stQTL 3aQTL enQTL
        do
                mkdir -p ${qtl}/${tis}/log ${qtl}/${tis}/summary ${qtl}/${tis}/covariances ${qtl}/${tis}/weights
		for chr in {1..26}
		do
                	jsub -q normal -n 1 -R "span[hosts=1]" -J train_${qtl}_${tis}_${chr} \
						-e ${qtl}/${tis}/log/03.train.${chr}.%J.log -o ${qtl}/${tis}/log/03.train.${chr}.%J.log \
						"Rscript 03.gtex_tiss_chrom_training.R ${chr} ${qtl}/${tis}"
		done
        done
done
rm -rf 3aQTL/Embryo
awk 'NF!=24{print FILENAME}' */*/summary/Model_training_chr*_model_summaries.txt > 03.error.summaries.list

# 4. Combine results
for tis in `cut -f1 tissue40.list | sed '1d'`
do
        for qtl in eQTL sQTL eeQTL isoQTL stQTL 3aQTL enQTL
        do
		mkdir -p ${qtl}/${tis}/dbs
                jsub -q normal -n 1 -R "span[hosts=1]" -J combine_${qtl}_${tis} \
					-e ${qtl}/${tis}/log/04.combine.%J.log -o ${qtl}/${tis}/log/04.combine.%J.log \
					"bash 04.combine.sh ${qtl} ${tis}"
        done
done
rm -rf 3aQTL/Embryo

# 5. generate parquet genotypes
mkdir -p parquet_genotypes/log
for tis in `cut -f1 tissue40.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J parquet_genotypes_${tis} \
		-e parquet_genotypes/log/05.parquet_genotypes.${tis}.%J.log -o parquet_genotypes/log/05.parquet_genotypes.${tis}.%J.log \
		"bash 05.parquet_genotypes.sh ${tis}"
done
## merge multi-tissue genotype files
mkdir -p parquet_genotypes/raw
ln -s /storage/public/home/2020060185/software/summary-gwas-imputation-master/src/genomic_tools_lib
for chr in {1..26}
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J variable_to_parquet_${chr} \
		-e parquet_genotypes/log/05.variable_to_parquet.chr${chr}.%J.log -o parquet_genotypes/log/05.variable_to_parquet.chr${chr}.%J.log \
		"bash 05.variable_to_parquet.sh ${chr}"
done

# 6. Construction of multi-tissue covariate files
for qtl in eQTL sQTL eeQTL isoQTL stQTL 3aQTL enQTL
do
	mkdir -p ${qtl}/log
	jsub -q normal -n 1 -R "span[hosts=1]" -J meta_${qtl} \
		-e ${qtl}/log/06.meta.%J.log -o ${qtl}/log/06.meta.%J.log \
		"bash 06.meta_covariance.sh ${qtl}"
done
