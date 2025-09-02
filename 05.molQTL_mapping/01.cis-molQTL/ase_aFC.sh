#!/bin/bash
tis=$1
threads=$2

mkdir -p ${tis}/results/phaser
# 1. preapare gw_phased.bed
echo "prepare gw_phased bed..."
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/03.MP2/ase/phaser_expr_matrix.gw_phased.bed.gz | cut -f1-4 > ${tis}/results/phaser/gw_phased.ann.bed
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/02.stat/03.MP2/ase/phaser_expr_matrix.gw_phased.bed.gz | cut -f5- | csvtk transpose -t -j $threads | csvtk join -t -j $threads -H -f '2;1' ${tis}/${tis}.list - | csvtk transpose -t -j $threads | sed '2d' > ${tis}/results/phaser/gw_phased.expr.bed
paste ${tis}/results/phaser/gw_phased.ann.bed ${tis}/results/phaser/gw_phased.expr.bed | bgzip -cf --threads $threads > ${tis}/results/phaser/${tis}.gw_phased.bed.gz

# 2. preapare vcf
echo "prepare vcf..."
cat ${tis}/${tis}.list | while read ind sample
do
    ls /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/08.ASE/${sample}/${sample}.vcf.gz
done > ${tis}/results/phaser/vcf.list
bcftools merge -l ${tis}/results/phaser/vcf.list -Oz -o ${tis}/results/phaser/${tis}.raw.vcf.gz --threads $threads
bcftools reheader -s ${tis}/${tis}.samplelist -o ${tis}/results/phaser/${tis}.phaser.vcf.gz ${tis}/results/phaser/${tis}.raw.vcf.gz --threads $threads
rm -f ${tis}/results/phaser/${tis}.raw.vcf.gz
bcftools index -t -f ${tis}/results/phaser/${tis}.phaser.vcf.gz --threads $threads

# 3. preapare pair and map file
echo "prepare map file ..."
paste ${tis}/${tis}.samplelist ${tis}/${tis}.samplelist | sed '1ivcf_sample\tbed_sample' > ${tis}/results/phaser/${tis}.map.txt
echo "prepare pair file ..."
## select variants with at least 10 individuals with ASE data and a minimum of 8 reads per individual
cat ${tis}/${tis}.list | while read ind sample
do
    sed '1d' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/02.MP2/01.pipeline/08.ASE/${sample}/${sample}.allelic_counts.txt | awk '$8>=8{print $1"\t"$2}'
done | sort -V | uniq -c | awk 'BEGIN{print "var_contig\tvar_pos"}$1>=10{print $2"\t"$3}' > ${tis}/results/phaser/allelic_counts.list
## select the top variant of each gene
zcat ${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | cut -f1,7 | csvtk join -t -H -f '2;2' - ${tis}/genotypes/${tis}.bim | cut -f1-3,5-7 | sed '1igene_id\tvar_id\tvar_contig\tvar_pos\tvar_alt\tvar_ref' > ${tis}/results/phaser/${tis}.pair.txt

# 4. run ase afc
source /storage/public/home/2020060185/anaconda3/envs/phaser/bin/activate phaser
python /storage/public/home/2020060185/software/phaser/phaser_pop/phaser_cis_var.py \
    --bed ${tis}/results/phaser/${tis}.gw_phased.bed.gz \
    --vcf ${tis}/results/phaser/${tis}.phaser.vcf.gz \
    --pair ${tis}/results/phaser/${tis}.pair.txt \
    --map ${tis}/results/phaser/${tis}.map.txt \
    --o ${tis}/results/phaser/${tis}.aFC.txt --t $threads
