#!/bin/bash
group1=Europe
group2=Central_and_East_Asia

# 1. get fd-eQTL
awk -v g=$group1 '$2==g&&$NF=="TRUE"{print $1"\t"$3}' tensorqtl_permutation.txt > group1.pairlist
awk -v g=$group2 '$2==g&&$NF=="TRUE"{print $1"\t"$3}' tensorqtl_permutation.txt > group2.pairlist
bash fd_eqtl.sh $group1 $group2

# 2. get ld-eQTL
## extract egenes
sort -u group1.pairlist group2.pairlist > egene.pairlist
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	mkdir -p ${tis}/phenotypes/mesusie
	grep -w "^${tis}" egene.pairlist | cut -f2 | split -dl 110
	mv x?? ${tis}/phenotypes/mesusie ; echo $tis
done
## run MESuSiE
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
    for num in `ls ${tis}/phenotypes/mesusie/x?? | sed 's/\//\t/g' | cut -f4`
    do
	jsub -q normal -n 3 -R "span[hosts=1]" -J v1.anc_mesusie_${tis}_${num} \
		-e ${tis}/log/05.mesusie.${tis}.${num}.%J.log -o ${tis}/log/05.mesusie.${tis}.${num}.%J.log \
		"bash mesusie.sh ${tis} ${num} ${group1} ${group2}"
    done
done
## combine results
for tis in `cut -f1 tissue.filtered2.list | sed '1d'`
do
	jsub -q normal -n 1 -R "span[hosts=1]" -J v1.anc_mesusie_combine_${tis} \
		-e ${tis}/log/05.mesusie_combine.${tis}.%J.log -o ${tis}/log/05.mesusie_combine.${tis}.%J.log \
		"bash mesusie_combine.sh $tis ${group1} ${group2}"
done
## combine all results
zcat */results/mesusie/*.mesusie.credible.group1.sig.gz | awk 'NR==1||$1!="tissue"' > mesusie.credible.group1.sig.txt
zcat */results/mesusie/*.mesusie.credible.group2.sig.gz | awk 'NR==1||$1!="tissue"' > mesusie.credible.group2.sig.txt
zcat */results/mesusie/*.mesusie.credible.shared.sig.gz | awk 'NR==1||$1!="tissue"' | gzip -c > mesusie.credible.shared.sig.gz
zcat */results/mesusie/*.mesusie.credible.shared.sig.lead.gz | awk 'NR==1||$1!="tissue"' | gzip -c > mesusie.credible.shared.sig.lead.gz

## get ld-eQTL
bash ld_eqtl.sh ${group1} ${group2}

# 3. get he-eQTL
bash he_eqtl.sh
