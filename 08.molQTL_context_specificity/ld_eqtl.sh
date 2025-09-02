#!/bin/bash
group1=$1
group2=$2

awk -v g=$group1 'NR==1||($1==g&&$11>=0.01&&$NF>=0.01)' susier.uniq.maf.txt > susier.group1.uniq.maf0.01.txt
awk -v g=$group2 'NR==1||($1==g&&$11>=0.01&&$NF>=0.01)' susier.uniq.maf.txt > susier.group2.uniq.maf0.01.txt

# 5. get ld-eQTLs (anc-specific PIP > 0.5)
awk '$7>0.5' mesusie.credible.group1.sig.txt | cut -f1-9 | csvtk join -t -f 'tissue,pheno_id,SNP;tissue,pheno_id,SNP' susier.group1.uniq.maf0.01.txt - | sed '1s/cs/cs_mesusie/2' > ld_eqtl.txt
awk '$8>0.5' mesusie.credible.group2.sig.txt | cut -f1-9 | csvtk join -t -f 'tissue,pheno_id,SNP;tissue,pheno_id,SNP' susier.group2.uniq.maf0.01.txt - | sed '1d' >> ld_eqtl.txt
awk '{if (NR==1){print}else{key=$1" "$2" "$3" "$5;if ($6>max[key]){max[key]=$6;line[key]=$0}}}END{for (k in line){print line[k]}}' ld_eqtl.txt > ld_eqtl.lead.txt
