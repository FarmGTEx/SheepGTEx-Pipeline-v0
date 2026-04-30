#!/bin/bash
for i in {1..3}
do
	mkdir -p lfsr${i}/tissues lfsr${i}/log
	for j in {3..53}
	do
		tis=`cut -f $j /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/00.lead/lfsr_sig_gene.count.txt | head -1`
		echo $i $tis
		awk -v i=${i} -v j=${j} '$2<=i&&$j=="True"{print $1}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/00.lead/lfsr_sig_gene.count.txt | sort -u | sed '1iphenotype_id' > lfsr${i}/tissues/${tis}.genelist
		if [[ -s "lfsr${i}/tissues/${tis}.genelist" ]]; then
			zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | csvtk join -t -f 'phenotype_id;phenotype_id' - lfsr${i}/tissues/${tis}.genelist | cut -f2 | sed '1d' | sort -u > lfsr${i}/tissues/${tis}.txt
  		else
    			echo "output file lfsr${i}/tissues/${tis}.genelist is empty!"
  		fi
		find lfsr${i}/tissues/${tis}.txt -type f -empty -delete
	done
done

for i in {49..51}
do
	mkdir -p lfsr${i}/tissues lfsr${i}/log
	for j in {3..53}
	do
		tis=`cut -f $j /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/00.lead/lfsr_sig_gene.count.txt | head -1`
		echo $i $tis
		awk -v i=${i} -v j=${j} '$2>=i&&$j=="True"{print $1}' /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/00.lead/lfsr_sig_gene.count.txt | sort -u | sed '1iphenotype_id' > lfsr${i}/tissues/${tis}.genelist
		if [[ -s "lfsr${i}/tissues/${tis}.genelist" ]]; then
			zcat /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v1.min40_split/${tis}/results/tensorqtl/nominal/${tis}.cis_qtl_pairs.sig.txt.gz | csvtk join -t -f 'phenotype_id;phenotype_id' - lfsr${i}/tissues/${tis}.genelist | cut -f2 | sed '1d' | sort -u > lfsr${i}/tissues/${tis}.txt
  		else
    			echo "output file lfsr${i}/tissues/${tis}.genelist is empty!"
  		fi
		find lfsr${i}/tissues/${tis}.txt -type f -empty -delete
	done
done

