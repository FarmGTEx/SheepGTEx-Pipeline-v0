#!/bin/bash
tis=$1

for qtl in sQTL eeQTL isoQTL stQTL 3aQTL enQTL
do
    echo $tis $qtl
    sed "s/^/${qtl}\t/g" ${tis}/results/coloc/eQTL_${qtl}.pph4 | sed "1s/^$qtl/qtl/" > ${tis}/results/coloc/eQTL_${qtl}.addcol.pph4
    zcat ${tis}/qtls/${qtl}.nominal/${tis}.cis_qtl_pairs.sig.txt.gz | cut -f1 | sed '1d' | sort -u > ${tis}/qtls/${qtl}.siglist
    if [ -s "${tis}/qtls/${qtl}.siglist" ]; then
        sed '1iphenotype2' ${tis}/qtls/${qtl}.siglist | csvtk join -t -f 'phenotype2;phenotype2' ${tis}/results/coloc/eQTL_${qtl}.addcol.pph4 - > ${tis}/results/coloc/eQTL_${qtl}.addcol.sig.pph4
    else
        echo "${tis}/qtls/${qtl}.siglist is empty!"
	rm -f ${tis}/results/coloc/eQTL_${qtl}.addcol.sig.pph4
    fi
done
awk 'NR==1||FNR>1' ${tis}/results/coloc/eQTL_*.addcol.pph4 > ${tis}/results/coloc/eQTL.addcol.pph4
awk 'NR==1||FNR>1' ${tis}/results/coloc/eQTL_*.addcol.sig.pph4 > ${tis}/results/coloc/eQTL.addcol.sig.pph4
