#!/bin/bash
tis=$1

mkdir -p ${tis}/qtls
# prepare permutations
## eQTL
zcat ../01.eQTL/v1.min40_split/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | awk '{print $1"\t"$0}' | sed '1s/phenotype_id/phenotype_id0/' > ${tis}/qtls/eQTL.permutation.txt
ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/nominal ${tis}/qtls/eQTL.nominal
## sQTL
zcat ../02.sQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | sed -e '1s/phenotype_id/phenotype_id0/' -e '1s/group_id/phenotype_id/' > ${tis}/qtls/sQTL.permutation.txt
ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/02.sQTL/${tis}/results/tensorqtl/nominal ${tis}/qtls/sQTL.nominal
## eeQTL
zcat ../03.eeQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | sed -e '1s/phenotype_id/phenotype_id0/' -e '1s/group_id/phenotype_id/' > ${tis}/qtls/eeQTL.permutation.txt
ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/${tis}/results/tensorqtl/nominal ${tis}/qtls/eeQTL.nominal
## isoQTL
zcat ../04.isoQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | sed -e '1s/phenotype_id/phenotype_id0/' -e '1s/group_id/phenotype_id/' > ${tis}/qtls/isoQTL.permutation.txt
ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/${tis}/results/tensorqtl/nominal ${tis}/qtls/isoQTL.nominal
## stQTL
zcat ../05.stQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | awk '{print $1"\t"$0}' | sed '1s/phenotype_id/phenotype_id0/' > ${tis}/qtls/stQTL.permutation.txt
ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/05.stQTL/${tis}/results/tensorqtl/nominal ${tis}/qtls/stQTL.nominal
## 3aQTL
zcat ../06.3aQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | sed -e '1s/phenotype_id/phenotype_id0/' -e '1s/group_id/phenotype_id/' > ${tis}/qtls/3aQTL.permutation.txt
ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/results/tensorqtl/nominal ${tis}/qtls/3aQTL.nominal
## enQTL
zcat ../07.enQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.txt.gz | sed -e '1s/phenotype_id/phenotype_id0/' -e '1s/group_id/phenotype_id/' > ${tis}/qtls/enQTL.permutation.txt
ln -s /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL/${tis}/results/tensorqtl/nominal ${tis}/qtls/enQTL.nominal
