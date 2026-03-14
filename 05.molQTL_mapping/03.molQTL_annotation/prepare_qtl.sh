#!/bin/bash
tis=$1

mkdir -p ${tis}/lead ${tis}/susier_lead
echo eQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/eQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/susier/${tis}.susier.credible.sig.lead.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier_lead/eQTL.txt

echo sQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/02.sQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/sQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/02.sQTL/${tis}/results/susier/${tis}.susier.credible.sig.lead.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier_lead/sQTL.txt

echo eeQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/eeQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/${tis}/results/susier/${tis}.susier.credible.sig.lead.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier_lead/eeQTL.txt

echo isoQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/isoQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/${tis}/results/susier/${tis}.susier.credible.sig.lead.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier_lead/isoQTL.txt

echo stQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/05.stQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/stQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/05.stQTL/${tis}/results/susier/${tis}.susier.credible.sig.lead.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier_lead/stQTL.txt

echo 3aQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/3aQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/results/susier/${tis}.susier.credible.sig.lead.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier_lead/3aQTL.txt

echo enQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.eenhancers.txt | sed '1d' | sort -u > ${tis}/lead/enQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL/${tis}/results/susier/${tis}.susier.credible.sig.lead.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier_lead/enQTL.txt

find ${tis}/lead -type f -empty -delete
find ${tis}/susier_lead -type f -empty -delete
