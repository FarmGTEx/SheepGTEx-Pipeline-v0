#!/bin/bash
## Effect size of eQTL
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J aFC0_${tis} \
        -e ${tis}/log/06.aFC0.${tis}.%J.log -o ${tis}/log/06.aFC0.${tis}.%J.log "bash aFC0.sh ${tis}"
done
# combine reasults
head -1 Heart/results/tensorqtl/permutation/Heart.log2aFC0.txt | sed 's/^/Tissue\t/g' > log2aFC0.txt
for tis in `cut -f1 ../tissue40.list | sed '1d'` ; do awk -v tis=$tis '{if (NR!=1) print tis"\t"$0}' ${tis}/results/tensorqtl/permutation/${tis}.log2aFC0.txt ; done >> log2aFC0.txt

# Measure cis-eQTL effect size by ASE
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        jsub -q normal -n 4 -R "span[hosts=1]" -J ASE_aFC_${tis} \
                -e ${tis}/log/06.ase_aFC.${tis}.%J.log -o ${tis}/log/06.ase_aFC.${tis}.%J.log "bash ase_aFC.sh ${tis} 4"
done
head -1 Heart/results/phaser/Heart.aFC.txt | csvtk cut -t -f 'gene,var_het_afc,var_het_pval' | sed 's/^/Tissue\t/g' > ase_aFC.txt
for tis in `cut -f1 ../tissue40.list | sed '1d'` ; do csvtk cut -t -f 'gene,var_het_afc,var_het_pval' ${tis}/results/phaser/${tis}.aFC.txt | awk -v tis=$tis '{if (NR!=1) print tis"\t"$0}' ; done >> ase_aFC.txt
