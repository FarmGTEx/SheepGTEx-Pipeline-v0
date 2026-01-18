chr=$1
ne=$2  # ne=1000
mkdir -p ./imputedbcf/${2}/

while IFS="" read -r LINE || [ -n "$LINE" ]; 
do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2)	
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
	/storage/public/home/2020060185/software/GLIMPSE2/GLIMPSE2_phase_static \
		--bam-list bamlist \
		--ne ${2} \
		--log ./logs/${2}.chr${CHR}_${REGS}_${REGE}.log \
		--reference /storage/public/home/2021050411/20.phase/sheep3125_snppanel/bin_refpanel/${chr}/Sheep3125.chr${chr}_${CHR}_${REGS}_${REGE}.bin \
		--output ./imputedbcf/${2}/chr${CHR}_${REGS}_${REGE}.bcf \
		--threads 30
done < /storage/public/home/2021050411/20.phase/sheep3125_snppanel/bin_refpanel/${chr}/chunks_chr${chr}.txt
chr=$1

mkdir -p ./imputedvcf

while IFS="" read -r LINE || [ -n "$LINE" ]; 
do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2)
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
	/storage/public/home/2020060185/software/GLIMPSE2/GLIMPSE2_ligate_static \
		--input ./imputedbcf/${chr}.ne${2}.bcflist \
		--output ./imputedvcf/${chr}.${2}.ligated.bcf \
		--threads 30
done < /storage/public/home/2021050411/20.phase/sheep3125_snppanel/bin_refpanel/${chr}/chunks_chr${chr}.txt


bcftools view ./imputedvcf/${chr}.${2}.ligated.bcf -Oz -o ./imputedvcf/${chr}.${2}.ligated.vcf.gz --threads 30
bcftools index -t ./imputedvcf/${chr}.${2}.ligated.vcf.gz --threads 30
chr=$1

mkdir -p recallINFO
source /storage/public/home/2020060185/.bashrc
#imputed vcf recall INFO
python recalINFO.py --invcf ./imputedvcf/${chr}.${2}.ligated.bcf | \
	bcftools +fill-tags -Oz -o ./recallINFO/chr${chr}.${2}.recalINFO.vcf.gz --threads 2 -- -t MAF
bcftools view -m2 -M2 ./recallINFO/chr${chr}.${2}.recalINFO.vcf.gz -Oz -o ./recallINFO/chr${chr}.${2}.snp.recalINFO.vcf.gz --threads 2
bcftools index -t ./recallINFO/chr${chr}.${2}.snp.recalINFO.vcf.gz --threads 2
bcftools filter -i 'MAF >= 0.01' ./recallINFO/chr${chr}.${2}.snp.recalINFO.vcf.gz -Oz -o ./recallINFO/chr${chr}.${2}.maf0.01.recalINFO.vcf.gz --threads 2
bcftools index -t ./recallINFO/chr${chr}.${2}.maf0.01.recalINFO.vcf.gz --threads 2

chr=$1
bcftools view -i 'INFO>=0.8' ./recallINFO/chr${chr}.${2}.maf0.01.recalINFO.vcf.gz --threads 4 -Oz -o info0.75/chr${chr}.${2}.maf0.01.filter.vcf.gz 
bcftools index info0.75/chr${chr}.${2}.maf0.01.filter.vcf.gz --threads 4
less -S chr${1}.1000.maf0.01.filter.vcf.gz |  awk 'BEGIN{OFS="\t"} /^#/ {print $0; next} {$3=$1"_"$2; print $0}' | bgzip -c > final/chr${1}.id.maf0.01.filter.vcf.gz
bcftools index final/chr${1}.id.maf0.01.filter.vcf.gz
bcftools concat \
	chr1.id.maf0.01.filter.vcf.gz \
	chr2.id.maf0.01.filter.vcf.gz \
	chr3.id.maf0.01.filter.vcf.gz \
	chr4.id.maf0.01.filter.vcf.gz \
	chr5.id.maf0.01.filter.vcf.gz \
	chr6.id.maf0.01.filter.vcf.gz \
	chr7.id.maf0.01.filter.vcf.gz \
	chr8.id.maf0.01.filter.vcf.gz \
	chr9.id.maf0.01.filter.vcf.gz \
	chr10.id.maf0.01.filter.vcf.gz \
	chr11.id.maf0.01.filter.vcf.gz \
	chr12.id.maf0.01.filter.vcf.gz \
	chr13.id.maf0.01.filter.vcf.gz \
	chr14.id.maf0.01.filter.vcf.gz \
	chr15.id.maf0.01.filter.vcf.gz \
	chr16.id.maf0.01.filter.vcf.gz \
	chr17.id.maf0.01.filter.vcf.gz \
	chr18.id.maf0.01.filter.vcf.gz \
	chr19.id.maf0.01.filter.vcf.gz \
	chr20.id.maf0.01.filter.vcf.gz \
	chr21.id.maf0.01.filter.vcf.gz \
	chr22.id.maf0.01.filter.vcf.gz \
	chr23.id.maf0.01.filter.vcf.gz \
	chr24.id.maf0.01.filter.vcf.gz \
	chr25.id.maf0.01.filter.vcf.gz \
	chr26.id.maf0.01.filter.vcf.gz \
	-Oz -o info0.75.maf0.01.filter.vcf.gz --threads 24

#predict gene tmm value
/storage/public/home/2022060212/anaconda3/envs/imlabtools/bin/python3 /storage/public/home/2022060212/software/MetaXcan-master/software/Predict.py \
--model_db_path /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/04.gwas_coloc/twas/01.makeModels/predictDB_elastic_net/eQTL/${1}/dbs/sheepgtex_models_filtered_signif.db \
--vcf_genotypes /storage/public/home/2021050411/20.phase/v3/info0.75/final/info0.75.maf0.01.filter.vcf.gz \
--vcf_mode genotyped \
--prediction_output out/0.5X.sum.predict.${1}.txt \
--prediction_summary_output out/0.5X.sum.${1}.txt \
--verbosity 9 \
--throw
#predict gene site frequency
bcftools view -R ${1}/${1}.ahcy.pos.skin final/info0.75.maf0.01.filter.vcf.gz -Oz -o ${1}/${1}.vcf.gz --threads 20
bcftools index ${1}/${1}.vcf.gz --threads 20
vcftools --gzvcf ${1}/${1}.vcf.gz --keep group/${2}.final.txt --freq --out ${1}/${1}.${2}.freq
