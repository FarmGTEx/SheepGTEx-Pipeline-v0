#!/bin/bash
group1=$1
group2=$2

# 1. find overlapped credible sets in each group
zcat susier.credible.sig.gz | awk -v g=$group1 'NR==1||$1==g' > susier.group1.txt
zcat susier.credible.sig.gz | awk -v g=$group2 'NR==1||$1==g' > susier.group2.txt
csvtk join -t -f 'tissue,pheno_id,SNP;tissue,pheno_id,SNP' susier.group1.txt susier.group2.txt > susier.overlap.txt
awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$18"\t"$2"\t"$3"\t"$19}' susier.overlap.txt | csvtk uniq -t -f '1,2,3,4,5,6,7,8' > susier.overlap.cs.list

# 2. find unique credible sets in each group
cut -f1-3,5 susier.group1.txt | csvtk uniq -t -f '1,2,3,4' > susier.group1.cs.list
cut -f1-3,5 susier.group2.txt | csvtk uniq -t -f '1,2,3,4' > susier.group2.cs.list
cut -f1-4 susier.overlap.cs.list | sort - susier.group1.cs.list | uniq -u | sed '1igroup\ttissue\tpheno_id\tcs' | csvtk join -t -f 'group,tissue,pheno_id,cs;group,tissue,pheno_id,cs' susier.group1.txt - > susier.group1.uniq.txt
cut -f5-8 susier.overlap.cs.list | sort - susier.group2.cs.list | uniq -u | sed '1igroup\ttissue\tpheno_id\tcs' | csvtk join -t -f 'group,tissue,pheno_id,cs;group,tissue,pheno_id,cs' susier.group2.txt - > susier.group2.uniq.txt

# 3. find MAF in another group
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	echo $tis
	# group1
	awk -v t=$tis 'NR==1||$2==t' susier.group1.uniq.txt > ${tis}/results/susier/${group1}/susier.uniq.txt
	zcat ${tis}/genotypes/${tis}.${group2}.frq.gz | awk '{print $2"\t"$5}' | csvtk join -t -f 'SNP;SNP' ${tis}/results/susier/${group1}/susier.uniq.txt - > ${tis}/results/susier/${group1}/susier.uniq.maf.txt
	## group2
	awk -v t=$tis 'NR==1||$2==t' susier.group2.uniq.txt > ${tis}/results/susier/${group2}/susier.uniq.txt
	zcat ${tis}/genotypes/${tis}.${group1}.frq.gz | awk '{print $2"\t"$5}' | csvtk join -t -f 'SNP;SNP' ${tis}/results/susier/${group2}/susier.uniq.txt - > ${tis}/results/susier/${group2}/susier.uniq.maf.txt
done

# 4. get fd-eQTLs
cat */results/susier/*/susier.uniq.maf.txt | awk 'NR==1||$1!="group"' > susier.uniq.maf.txt
awk 'NR==1||($11>=0.05&&$NF<0.01)' susier.uniq.maf.txt | sed -i '1s/MAF/MAF.1/2' > fd_eqtl.txt
awk '{if (NR==1){print}else{key=$1" "$2" "$3" "$5;if ($6>max[key]){max[key]=$6;line[key]=$0}}}END{for (k in line){print line[k]}}' fd_eqtl.txt > fd_eqtl.lead.txt
