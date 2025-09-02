#!/bin/bash
tis=$1
group=$2
num=$3

mkdir -p ${tis}/results/susier/${group}/${num}
# prepare summary data
for gene in `cat ${tis}/phenotypes/${group}/${num}`
do
  mkdir -p ${tis}/local/susier/${group}/${gene}
  chr=`zcat ${tis}/phenotypes/${tis}.${group}.expression.bed.gz | awk -v g=${gene} '$4==g' | cut -f1`
  sample_size=`cat ${tis}/genotypes/${tis}.${group}.fam | wc -l`
  # prepare summary statistics
  zcat ${tis}/results/tensorqtl/nominal/${tis}.${group}.cis_qtl_pairs.${chr}.txt.gz | awk -v g=${gene} '$1==g&&NF==9' | csvtk join -t -H -f '2;2' ${tis}/genotypes/${tis}.${group}.bim - | awk -v OFS="\t" -v n=$sample_size 'BEGIN{print "CHR_POS\tCHR\tPOS\tREF\tALT\tMAF\tBETA\tSE\tZ\tN"}{if ($9 <= 0.5) print $2,$1,$4,$6,$5,$9,$13,$14,$13/$14,n; else print $2,$1,$4,$6,$5,1-$9,$13,$14,$13/$14,n}' > ${tis}/local/susier/${group}/${gene}/summary.txt
  # prepare genotype
  cut -f1 ${tis}/local/susier/${group}/${gene}/summary.txt | sed '1d' > ${tis}/local/susier/${group}/${gene}/snp.list
  plink --bfile ${tis}/genotypes/${tis}.${group} --extract ${tis}/local/susier/${group}/${gene}/snp.list --sheep --keep-allele-order --make-bed --out ${tis}/local/susier/${group}/${gene}/genotype
  # SuSiE
  Rscript /storage/public/home/2020060185/script/SuSiE.R \
    ${tis}/local/susier/${group}/${gene}/summary.txt \
    ${tis}/local/susier/${group}/${gene}/genotype \
    ${tis}/local/susier/${group}/${gene}/${tis}.${gene}
  date
  # save variants
  if [[ -s "${tis}/local/susier/${group}/${gene}/${tis}.${gene}.txt" ]]; then
    awk -v t=$tis -v g=$gene '{if (NR==1){print "tissue\tpheno_id\t"$0}else{print t"\t"g"\t"$0}}' ${tis}/local/susier/${group}/${gene}/${tis}.${gene}.txt | csvtk join -t -f 'SNP;CHR_POS' - ${tis}/local/susier/${group}/${gene}/summary.txt | gzip -c > ${tis}/results/susier/${group}/${num}/${tis}.${gene}.susier.gz
  else
    echo "output file ${tis}/local/susier/${group}/${gene}/${tis}.${gene}.txt does not exists!"
  fi
  rm -rf ${tis}/local/susier/${group}/${gene}
done

# combine results
zcat ${tis}/results/susier/${group}/${num}/${tis}.*.susier.gz | awk 'NR==1||$1!="tissue"' | gzip -c > ${tis}/results/susier/${group}/${num}.susier.gz
rm -rf ${tis}/results/susier/${group}/${num}
