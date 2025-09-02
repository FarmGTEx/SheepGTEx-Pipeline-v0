#!/bin/bash
tis=$1
num=$2
threads=$3

mkdir -p ${tis}/opera/output/${num}
for gene in `cat ${tis}/opera/split/${num}`
do
	mkdir -p ${tis}/opera/local/${gene}
	chr=`zcat ${tis}/phenotypes/eQTL.expression.bed.gz | awk -v g=${gene} '$4==g' | cut -f1`
	pthresh=`awk -v g=${gene} '$1==g{print $(NF-1)}' ${tis}/qtls/eQTL.permutation.txt`
	lead=`awk -v g=${gene} '$1==g{print $7}' ${tis}/qtls/eQTL.permutation.txt`

	# stage 0
	## prepare eQTL input as outcome in GCTA-COJO format
	zcat ${tis}/qtls/eQTL.nominal/${tis}.cis_qtl_pairs.${chr}.txt.gz | awk -v g=${gene} 'NF==9&&$1==g' | awk -v OFS="\t" '{if($4<0.5){print $2,$4,$8,$9,$7,$6/$4/2}else{print $2,$4,$8,$9,$7,$6/(1-$4)/2}}' | csvtk join -t -H -j ${threads} -f '1;2' - ${tis}/opera/input/${tis}.bim | awk -v OFS="\t" '{print $1,$10,$11,$2,$3,$4,$5,$6}' | sed '1iSNP\tA1\tA2\tfrq\tb\tse\tP\tN' > ${tis}/opera/local/${gene}/eQTL.ma
	cut -f1 ${tis}/opera/local/${gene}/eQTL.ma | sed '1d' > ${tis}/opera/local/${gene}/cis.snplist
	## extract the LD independent eQTL loci from independent QTLs (add lead SNP if disappeared in some genes)
	zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/independent/${tis}.cis_independent_qtl.txt.gz | awk -v g=${gene} -v l=${lead} 'BEGIN{print l}$1==g{print $7}' | sort -u | sed 's/_/\t/g' | awk '{print$1"\t"$1"_"$2"\t"$2}' | sed '1iChr\tSNP\tbp' > ${tis}/opera/local/${gene}/cis.independent.loci
	site_num=`sed '1d' ${tis}/opera/local/${gene}/cis.independent.loci | wc -l`
       	echo "$gene eQTL done."
	
	# stage 1
	opera_Linux --besd-flist ${tis}/opera/input/besd.${chr}.list --gwas-summary ${tis}/opera/local/${gene}/eQTL.ma \
		--bfile ${tis}/opera/input/${tis} --sample-overlap --estimate-pi --out ${tis}/opera/local/${gene}/stage1 --thread-num $threads \
		--diff-freq 0.8 --diff-freq-prop 0.2 --peqtl-smr $pthresh --extract-snp ${tis}/opera/local/${gene}/cis.snplist
	echo "$gene stage1 done."

	if [ -s "${tis}/opera/local/${gene}/stage1.rho" ]; then	
		# stage 2
		opera_Linux --besd-flist ${tis}/opera/input/besd.${chr}.list --gwas-summary ${tis}/opera/local/${gene}/eQTL.ma \
			--bfile ${tis}/opera/input/${tis} --sample-overlap --rho-file ${tis}/opera/local/${gene}/stage1.rho --print-combo-ppa-res \
			--prior-pi-file ${tis}/opera/local/${gene}/stage1.pi --prior-var-file ${tis}/opera/local/${gene}/stage1.var \
			--out ${tis}/opera/local/${gene}/stage2 --thread-num $threads --extract-snp ${tis}/opera/local/${gene}/cis.snplist \
			--extract-gwas-loci ${tis}/opera/local/${gene}/cis.independent.loci
		echo "$gene stage2 done."

	 	# save results
		awk -v t=$tis -v g=$gene '{if(NR==1){print "Tissue\tGene\t"$0}else{print t"\t"g"\t"$0}}' ${tis}/opera/local/${gene}/stage2_combo.res | pigz -p $threads -c > ${tis}/opera/output/${num}/${gene}.stage2_combo.res.gz
		awk -v t=$tis -v g=$gene '{if(NR==1){print "Tissue\tGene\t"$0}else{print t"\t"g"\t"$0}}' ${tis}/opera/local/${gene}/stage2.prop | pigz -p $threads -c > ${tis}/opera/output/${num}/${gene}.stage2.prop.gz
		echo -e "$tis\t$gene\t$site_num" > ${tis}/opera/output/${num}/${gene}.sites
		for i in {1..6} ; do awk -v t=$tis -v g=$gene '{if(NR==1){print "Tissue\tGene\t"$0}else{print t"\t"g"\t"$0}}' ${tis}/opera/local/${gene}/stage2_${i}_expos_assoc.ppa | pigz -p $threads -c > ${tis}/opera/output/${num}/${gene}.stage2_${i}_expos_assoc.ppa.gz ; done
	else
		echo "$gene rho could not be calculated!"
	fi		
	rm -rf ${tis}/opera/local/${gene}
done

# combine results
zcat ${tis}/opera/output/${num}/*.stage2_combo.res.gz | awk 'NR==1||$1!="Tissue"' | pigz -p $threads -c > ${tis}/opera/output/${num}.stage2_combo.res.gz
zcat ${tis}/opera/output/${num}/*.stage2.prop.gz | awk 'NR==1||$1!="Tissue"' | pigz -p $threads -c > ${tis}/opera/output/${num}.stage2.prop.gz
cat ${tis}/opera/output/${num}/*.sites | pigz -p $threads -c > ${tis}/opera/output/${num}.sites.gz
for i in {1..6} ; do echo $i ; zcat ${tis}/opera/output/${num}/*.stage2_${i}_expos_assoc.ppa.gz | awk 'NR==1||$1!="Tissue"' | pigz -p $threads -c > ${tis}/opera/output/${num}.stage2_${i}_expos_assoc.ppa.gz ; done
rm -rf ${tis}/opera/output/${num}

