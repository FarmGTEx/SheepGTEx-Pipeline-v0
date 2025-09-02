#!/bin/bash
tis1=$1
tis2=$2

# combine results
sed '1d' combinations/${tis1}_${tis2}/egene.joined.permutation.txt | sed '1itissue1\tphenotype_id\tvariant_id1\tpval_nominal1\tslope1\tslope_se1\tis_eGene1\ttissue2\tvariant_id2\tpval_nominal2\tslope2\tslope_se2\tis_eGene2' | csvtk join -t -L -f 'phenotype_id;phenotype1' - combinations/${tis1}_${tis2}/coloc.pph4 | csvtk join -t -L -f 'phenotype_id;gene' - combinations/${tis1}_${tis2}/plink.r2 > combinations/${tis1}_${tis2}/combined.txt

# tissue specific eQTLs classification
python tissue_specific_effect_qtl.py \
	--inlfsr /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/00.lead/output/top_pairs/lfsr_m.s.RDS \
	--infile combinations/${tis1}_${tis2}/combined.txt \
	--outfile combinations/${tis1}_${tis2}/oppo_effect_eqtl.txt
rm -f combinations/${tis1}_${tis2}/combined.txt
