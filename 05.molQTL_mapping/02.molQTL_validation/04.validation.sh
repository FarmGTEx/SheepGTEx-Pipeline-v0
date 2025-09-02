#!/bin/bash
## invernal validation
mkdir -p log
for tis in `cut -f1 tissue100.list | sed '1d'`
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J internal_validation_${tis} \
                -e log/04.validation.${tis}.%J.log -o log/04.validation.${tis}.%J.log \
                "Rscript pi0est.R ${tis}/group1/results/tensorqtl ${tis}/group2/results/tensorqtl group1 group2 ${tis} ${tis}/${tis}"
done

# combine results
cat */*.pi1.txt | sed '1iDiscovery\tValidation\tTissue\tPi1\tSlope\tZ'> pi1.txt
