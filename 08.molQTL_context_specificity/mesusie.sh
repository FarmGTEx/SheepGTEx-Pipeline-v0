#!/bin/bash
tis=$1
num=$2
group1=$3
group2=$4

mkdir -p ${tis}/results/mesusie/${num}
# prepare summary data
for gene in `cat ${tis}/phenotypes/mesusie/${num}`
do
  mkdir -p ${tis}/local/mesusie/${gene}
  for group in ${group1} ${group2}
  do
    chr=`zcat ${tis}/phenotypes/${tis}.${group}.expression.bed.gz | awk -v g=${gene} '$4==g' | cut -f1`
    sample_size=`cat ${tis}/genotypes/${tis}.${group}.fam | wc -l`
    # prepare summary statistics
    zcat ${tis}/results/tensorqtl/nominal/${tis}.${group}.cis_qtl_pairs.${chr}.txt.gz | awk -v g=${gene} '$1==g&&NF==9' | csvtk join -t -H -f '2;2' ${tis}/genotypes/${tis}.${group}.bim - | awk -v OFS="\t" -v n=$sample_size 'BEGIN{print "CHR_POS\tCHR\tPOS\tREF\tALT\tMAF\tBETA\tSE\tZ\tN"}{if ($9 <= 0.5) print $2,$1,$4,$6,$5,$9,$13,$14,$13/$14,n; else print $2,$1,$4,$6,$5,1-$9,$13,$14,$13/$14,n}' > ${tis}/local/mesusie/${gene}/${group}.summary.txt
    # prepare genotype
    cut -f1 ${tis}/local/mesusie/${gene}/${group}.summary.txt | sed '1d' > ${tis}/local/mesusie/${gene}/${group}.snp.list
    plink --bfile ${tis}/genotypes/${tis}.${group} --extract ${tis}/local/mesusie/${gene}/${group}.snp.list --sheep --keep-allele-order --make-bed --out ${tis}/local/mesusie/${gene}/${group}
  done
  # MESuSiE
  Rscript MESuSiE.R \
    ${tis}/local/mesusie/${gene}/${group1}.summary.txt ${tis}/local/mesusie/${gene}/${group2}.summary.txt \
    ${tis}/local/mesusie/${gene}/${group1} ${tis}/local/mesusie/${gene}/${group2} \
    ${tis}/local/mesusie/${gene}/${tis}.${gene}
  date
  # save variants
  if [[ -s "${tis}/local/mesusie/${gene}/${tis}.${gene}.txt" ]]; then
    csvtk join -t -f 'CHR_POS,CHR,POS,REF,ALT;CHR_POS,CHR,POS,REF,ALT' ${tis}/local/mesusie/${gene}/${group1}.summary.txt ${tis}/local/mesusie/${gene}/${group2}.summary.txt > ${tis}/local/mesusie/${gene}/summary.txt
    awk -v t=$tis -v g=$gene '{if (NR==1){print "tissue\tpheno_id\t"$0}else{print t"\t"g"\t"$0}}' ${tis}/local/mesusie/${gene}/${tis}.${gene}.txt | csvtk join -t -f 'SNP;CHR_POS' - ${tis}/local/mesusie/${gene}/summary.txt | gzip -c > ${tis}/results/mesusie/${num}/${tis}.${gene}.mesusie.gz
    mv ${tis}/local/mesusie/${gene}/${tis}.${gene}.png ${tis}/results/mesusie/${num}
  else
    echo "output file ${tis}/local/mesusie/${gene}/${tis}.${gene}.txt does not exists!"
  fi
  rm -rf ${tis}/local/mesusie/${gene}
done

# combine results
zcat ${tis}/results/mesusie/${num}/${tis}.*.mesusie.gz | awk 'NR==1||$1!="tissue"' | gzip -c > ${tis}/results/mesusie/${num}.mesusie.gz
rm -f ${tis}/results/mesusie/${num}/${tis}.*.mesusie.gz
