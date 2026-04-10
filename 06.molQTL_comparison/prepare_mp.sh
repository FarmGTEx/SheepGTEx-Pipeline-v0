#!/bin/bash
tis=$1

mkdir -p ${tis}/phenotypes
# prepare mol-phenotype
## eQTL
cp ../01.eQTL/v1.min40_split/${tis}/phenotypes/${tis}.expression.bed.gz ${tis}/phenotypes/eQTL.expression.bed.gz
## sQTL
sed '1iphenotype_id\tgene_id' ../02.sQTL/${tis}/phenotypes/${tis}.leafcutter.phenotype_groups.txt > ${tis}/phenotypes/sQTL.phenotype_groups.txt
zcat ../02.sQTL/${tis}/phenotypes/${tis}.leafcutter.sorted.bed.gz | csvtk join -t -C '$' -f 'ID;phenotype_id' - ${tis}/phenotypes/sQTL.phenotype_groups.txt | sed '1s/ID/phenotype_id/' | bgzip -cf > ${tis}/phenotypes/sQTL.expression.bed.gz
## eeQTL
sed '1iphenotype_id\tgene_id' ../03.eeQTL/${tis}/phenotypes/${tis}.phenotype_groups.txt > ${tis}/phenotypes/eeQTL.phenotype_groups.txt
zcat ../03.eeQTL/${tis}/phenotypes/${tis}.expression.bed.gz | csvtk join -t -C '$' -f 'phenotype_id;phenotype_id' - ${tis}/phenotypes/eeQTL.phenotype_groups.txt | bgzip -cf > ${tis}/phenotypes/eeQTL.expression.bed.gz
## isoQTL
sed '1iphenotype_id\tgene_id' ../04.isoQTL/${tis}/phenotypes/${tis}.phenotype_groups.txt > ${tis}/phenotypes/isoQTL.phenotype_groups.txt
zcat ../04.isoQTL/${tis}/phenotypes/${tis}.expression.bed.gz | csvtk join -t -C '$' -f 'phenotype_id;phenotype_id' - ${tis}/phenotypes/isoQTL.phenotype_groups.txt | bgzip -cf > ${tis}/phenotypes/isoQTL.expression.bed.gz
## stQTL
cp ../05.stQTL/${tis}/phenotypes/${tis}.expression.bed.gz ${tis}/phenotypes/stQTL.expression.bed.gz
## 3aQTL
sed '1iphenotype_id\tgene_id' ../06.3aQTL/${tis}/phenotypes/${tis}.phenotype_groups.txt > ${tis}/phenotypes/3aQTL.phenotype_groups.txt
zcat ../06.3aQTL/${tis}/phenotypes/${tis}.expression.bed.gz | csvtk join -t -C '$' -f 'ID;phenotype_id' - ${tis}/phenotypes/3aQTL.phenotype_groups.txt | sed '1s/ID/phenotype_id/' | bgzip -cf > ${tis}/phenotypes/3aQTL.expression.bed.gz
## enQTL
sed '1iphenotype_id\tgene_id' ../07.enQTL_new/enQTL_new/${tis}/phenotypes/${tis}.phenotype_groups.txt > ${tis}/phenotypes/enQTL.phenotype_groups.txt
zcat ../07.enQTL_new/enQTL_new/${tis}/phenotypes/${tis}.all.expression.bed.gz | csvtk join -t -C '$' -f 'phenotype_id;phenotype_id' - ${tis}/phenotypes/enQTL.phenotype_groups.txt | bgzip -cf > ${tis}/phenotypes/enQTL.expression.bed.gz
