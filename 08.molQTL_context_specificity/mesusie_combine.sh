#!/bin/bash
tis=$1
group1=$2
group2=$3

# combine mesusie results
zcat ${tis}/results/mesusie/x??.mesusie.gz | awk 'NR==1||($1!="tissue"&&NF>=21)' | gzip -c > ${tis}/results/mesusie/${tis}.mesusie.gz
# extract 95% credible sets
zcat ${tis}/results/mesusie/${tis}.mesusie.gz | awk -v FS="\t" '$4!=""' | gzip -c > ${tis}/results/mesusie/${tis}.mesusie.credible.gz
# extract group1-specific eQTLs in 95% credible sets
zcat ${tis}/results/mesusie/${tis}.mesusie.credible.gz | awk 'NR==1||$5=="EUR"' | csvtk join -t -f 'pheno_id,SNP;phenotype_id,variant_id' - <(zcat ${tis}/results/tensorqtl/nominal/${tis}.${group1}.cis_qtl_pairs.sig.txt.gz | cut -f1,2,7,10) | gzip -c > ${tis}/results/mesusie/${tis}.mesusie.credible.group1.sig.gz
# extract group2-specific eQTLs in 95% credible sets
zcat ${tis}/results/mesusie/${tis}.mesusie.credible.gz | awk 'NR==1||$5=="CEA"' | csvtk join -t -f 'pheno_id,SNP;phenotype_id,variant_id' - <(zcat ${tis}/results/tensorqtl/nominal/${tis}.${group2}.cis_qtl_pairs.sig.txt.gz | cut -f1,2,7,10) | gzip -c > ${tis}/results/mesusie/${tis}.mesusie.credible.group2.sig.gz
# extract shared eQTLs (in at least one group) in 95% credible sets
zcat ${tis}/results/tensorqtl/nominal/${tis}.${group1}.cis_qtl_pairs.sig.txt.gz | cut -f1,2 | sed '1d' | sort -u > ${tis}/results/mesusie/group1.siglist
zcat ${tis}/results/tensorqtl/nominal/${tis}.${group2}.cis_qtl_pairs.sig.txt.gz | cut -f1,2 | sed '1d' | sort -u > ${tis}/results/mesusie/group2.siglist
sort -u ${tis}/results/mesusie/group1.siglist ${tis}/results/mesusie/group2.siglist | sed '1ipheno_id\tSNP' > ${tis}/results/mesusie/qtl.siglist
zcat ${tis}/results/mesusie/${tis}.mesusie.credible.gz | awk 'NR==1||$5=="EUR_CEA"' | csvtk join -t -f 'pheno_id,SNP;pheno_id,SNP' - ${tis}/results/mesusie/qtl.siglist | gzip -c > ${tis}/results/mesusie/${tis}.mesusie.credible.shared.sig.gz
# extract the lead eQTL in each shared 95% credible set
zcat ${tis}/results/mesusie/${tis}.mesusie.credible.shared.sig.gz | awk '{if (NR==1){print}else{key=$2" "$4;if ($6>max[key]){max[key]=$6;line[key]=$0}}}END{for (k in line){print line[k]}}' | gzip -c > ${tis}/results/mesusie/${group}/${tis}.mesusie.credible.shared.sig.lead.gz
