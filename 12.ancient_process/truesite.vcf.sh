###00.mplieup.sh 
chr=$1
bcftools mpileup \
	-O z -A \
   	-b bamlist \
        -q 30 \
        -Q 30 \
	-f ./reference/ARS-UI_Ramb_v2.0_CM2CHR.fasta \
	--annotate  FORMAT/DP,FORMAT/AD,INFO/AD \
	--threads 8 \
	-r ${chr} --ignore-RG  >  mplieup1/chr${chr}.mpilup.vcf.gz
bcftools index -t mplieup1/chr${1}.mpilup.vcf.gz --threads 8

####01.call_var_new.sh 
chr=$1

bcftools call \
-m \
-f GQ,GP \
-r ${chr} \
-T ./sheep3215.site.chr${chr}.vcf \
--threads 10 \
--skip-variants indels \
-Oz -o  call_vcf1/chr${chr}.vcf.gz mplieup1/chr${chr}.mpilup.vcf.gz
bcftools index -t call_vcf1/chr${chr}.vcf.gz --threads 10

###02.sample.sh 
bcftools view -s ${1} -Oz -o sample/${1}.chr${2}.vcf.gz out5/chr${2}.reheader.vcf.gz --threads 2
bcftools index sample/${1}.chr${2}.vcf.gz --threads 2

###03.filter.sh
bcftools filter \
  -S . \
  -e ' FMT/GQ<30 || FMT/DP<3 || FMT/DP>27 ' sample/${1}.chr${2}.vcf.gz \ ## if sample depth =9, choose DP >3 and DP<27
  -Oz -o filter/${1}.chr${2}.vcf.gz --threads 4
bcftools index filter/${1}.chr${2}.vcf.gz --threads 4

##5.reheader.sh
bcftools reheader -s rename.sample call_vcf1/chr${1}.vcf.gz -o out5/chr${1}.reheader.vcf.gz
bcftools index out5/chr${1}.reheader.vcf.gz


