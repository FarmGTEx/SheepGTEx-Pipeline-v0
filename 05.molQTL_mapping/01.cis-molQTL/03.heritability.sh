#!/bin/bash

# heritability (omiga)
for tis in `cut -f1 ../tissue40.list | sed '1d'`
do
        jsub -q normal -n 4 -R "span[hosts=1]" -J omiga_heritability_${tis} \
                -e ${tis}/log/03.omiga_heritability.${tis}.%J.log -o ${tis}/log/03.omiga_heritability.${tis}.%J.log \
                "bash omiga_heritability.sh ${tis} 4"
done

# combine results
zcat Heart/results/omiga/heritability/Heart.cis.her_est.txt.gz | head -1 | sed 's/^/Tissue\tSample_size\t/g' > omiga_h2.single.txt
sed '1d' ../tissue40.list | while read tis size num ; do zcat ${tis}/results/omiga/heritability/${tis}.cis.her_est.txt.gz | grep -v '^pheno_id' | awk -v tis=$tis -v size=$size '{print tis"\t"size"\t"$0}' ; zcat ${tis}/results/omiga/heritability/${tis}.trans.her_est.txt.gz | grep -v '^pheno_id' | awk -v tis=$tis -v size=$size '{print tis"\t"size"\t"$0}' ; done >> omiga_h2.single.txt
