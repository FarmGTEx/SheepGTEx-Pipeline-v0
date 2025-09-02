#!/bin/bash
tis=$1

mkdir -p ${tis}/lead ${tis}/qtls ${tis}/susier
echo eQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/eQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f2 | sed '1d' | sort -u > ${tis}/qtls/eQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/susier/${tis}.susier.credible.sig.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier/eQTL.txt

echo sQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/02.sQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/sQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/02.sQTL/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f2 | sed '1d' | sort -u > ${tis}/qtls/sQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/02.sQTL/${tis}/results/susier/${tis}.susier.credible.sig.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier/sQTL.txt

echo eeQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/eeQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f2 | sed '1d' | sort -u > ${tis}/qtls/eeQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/03.eeQTL/${tis}/results/susier/${tis}.susier.credible.sig.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier/eeQTL.txt

echo isoQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/isoQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f2 | sed '1d' | sort -u > ${tis}/qtls/isoQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/04.isoQTL/${tis}/results/susier/${tis}.susier.credible.sig.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier/isoQTL.txt

echo stQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/05.stQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/stQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/05.stQTL/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f2 | sed '1d' | sort -u > ${tis}/qtls/stQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/05.stQTL/${tis}/results/susier/${tis}.susier.credible.sig.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier/stQTL.txt

echo 3aQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.egenes.txt | sed '1d' | sort -u > ${tis}/lead/3aQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f2 | sed '1d' | sort -u > ${tis}/qtls/3aQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/06.3aQTL/${tis}/results/susier/${tis}.susier.credible.sig.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier/3aQTL.txt

echo enQTL
cut -f7 /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL_new/enQTL_new/${tis}/results/tensorqtl/permutation/${tis}.cis_qtl_fdr0.05.eenhancers.txt | sed '1d' | sort -u > ${tis}/lead/enQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL_new/enQTL_new/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f2 | sed '1d' | sort -u > ${tis}/qtls/enQTL.txt
zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/07.enQTL_new/enQTL_new/${tis}/results/susier/${tis}.susier.credible.sig.gz | cut -f3 | sed '1d' | sort -u > ${tis}/susier/enQTL.txt

# remove empty files
find ${tis}/lead -type f -empty -delete
find ${tis}/qtls -type f -empty -delete
find ${tis}/susier -type f -empty -delete
