### downsample.sh
samtools view -s ${1} ./${2}.filtered.bam -b -o ${2}.${3}X.nuclear.bam -@ 2

samtools index -b ${2}.${3}X.nuclear.bam -@ 2

#1.impute_bcf.sh 
chr=$1

mkdir -p ./imputedbcf1/${2}/

while IFS="" read -r LINE || [ -n "$LINE" ]; 
do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2)	
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
./GLIMPSE2_phase_static \
		--bam-list ${2}X.bamlist \
		--ne 1000 \
		--log ./log1/${2}.chr${CHR}_${REGS}_${REGE}.log \
		--reference ./sheep_panel/v1/bin_refpanel/${chr}/Sheep3125.chr${chr}_${CHR}_${REGS}_${REGE}.bin \
		--output ./imputedbcf1/${2}/chr${CHR}_${REGS}_${REGE}.bcf \
		--threads 25
done < ./bin_refpanel/${chr}/chunks_chr${chr}.txt

###2.bcf_vcf.sh 
chr=$1

mkdir -p ./imputedvcf1/${2}/

while IFS="" read -r LINE || [ -n "$LINE" ]; 
do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2)
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
./GLIMPSE2/GLIMPSE2_ligate_static \
		--input ./imputedbcf1/${chr}.${2}X.bcflist \
		--output ./imputedvcf1/${2}/${chr}.${2}X.ligated.bcf \
		--threads 30
done < ./sheep_panel/v1/bin_refpanel/${chr}/chunks_chr${chr}.txt


bcftools view ./imputedvcf1/${2}/${chr}.${2}X.ligated.bcf -Oz -o ./imputedvcf1/${2}/${chr}.${2}X.ligated.vcf.gz --threads 30
bcftools index -t ./imputedvcf1/${2}/${chr}.${2}X.ligated.vcf.gz --threads 30


###3.recalINFO.sh 
#!/bin/bash
chr=$1

mkdir -p recallINFO/${2}/
mkdir -p accuracy/${2}/
#imputed vcf recall INFO
python3 ./Script/recalINFO.py \
	--invcf ./imputedvcf1/${2}/${chr}.${2}X.ligated.bcf | \
bcftools +fill-tags -Oz -o ./recallINFO/${2}/chr${chr}.${2}X.recalINFO.vcf.gz --threads 2 -- -t MAF
bcftools view -m2 -M2 ./recallINFO/${2}/chr${chr}.${2}X.recalINFO.vcf.gz -Oz -o ./recallINFO/${2}/chr${chr}.${2}X.snp.recalINFO.vcf.gz --threads 2
bcftools index -t ./recallINFO/${2}/chr${chr}.${2}X.snp.recalINFO.vcf.gz --threads 2


##imputation accuracy 
#target samples rna imputation vcf
bcftools query -S samplelist -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF\t%INFO/MAF\t%INFO/INFO[\t%GT]\n' ./recallINFO/${2}/chr${chr}.${2}X.snp.recalINFO.vcf.gz | sed -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' -e 's/\.\/\.//g' -e 's/# //g' -e 's/:GT//g' -e 's/\[[^][]*\]//g' > accuracy/${2}/${chr}.${2}X.txt

###4.sample.sh 
chr=$1
bcftools query  -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF\t%INFO/MAF\t%INFO/INFO[\t%GT]\n' ./recallINFO/${2}/chr${chr}.${2}X.snp.recalINFO.vcf.gz | sed -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' -e 's/\.\/\.//g' -e 's/# //g' -e 's/:GT//g' -e 's/\[[^][]*\]//g' > accuracy/${2}/${chr}.${2}X.txt



####5.call.sh 
#cut -f1-3 accuracy/${1}_wgs0.txt | csvtk join -t -f 'CHROM,POS,REF;CHROM,POS,REF' - accuracy/${2}/${1}.${2}X.txt > accuracy1/${1}_rna.${2}X.txt
#cut -f1-3 accuracy/${2}/${1}.${2}X.txt | csvtk join -t -f 'CHROM,POS,REF;CHROM,POS,REF' - accuracy/${1}_wgs0.txt > accuracy1/${1}_wgs.${2}X.txt
# accuracy
python ./script/impute_accuracy.py \
	--imp_wgs list.${2}X \
	--wgsfile accuracy1/${1}_wgs.${2}X.txt \
	--impfile accuracy1/${1}_rna.${2}X.txt \
	--samplefile accuracy2/${1}.${2}X.sample.accuracy \
	--sitefile accuracy2/${1}.${2}X.sites.accuracy

##ts.sh 
bcftools view -i '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="G" & ALT="A") | (REF="A" & ALT="G")' sheep2.chr${1}.maf0.01miss0.1.vcf.gz --threads 4 -Oz -o sheep.ts.chr${1}.vcf.gz
bcftools index sheep.ts.chr${1}.vcf.gz --threads 4
##tv.sh 
bcftools index sheep2.chr${1}.maf0.01miss0.1.vcf.gz --threads 4
bcftools view \
	-e '(REF="C" & ALT="T") | (REF="T" & ALT="C") | (REF="G" & ALT="A") | (REF="A" & ALT="G")' sheep2.chr${1}.maf0.01miss0.1.vcf.gz --threads 4 \
       	-Oz -o sheep.tv.chr${1}.vcf.gz
bcftools index sheep.tv.chr${1}.vcf.gz --threads 4


