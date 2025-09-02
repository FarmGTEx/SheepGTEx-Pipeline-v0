#!/bin/bash
python he_qtl.py --susier_file susier.overlap.txt --mesusie_file mesusie.credible.shared.sig.gz --outfile he_eqtl.all.txt
awk 'NR==1||$NF=="True"' he_eqtl.all.txt > he_eqtl.txt
awk 'NR==1||($8>0&&$13<0)||($8<0&&$13>0)' he_eqtl.txt > he_eqtl.oppo.txt
## get lead eQTLs
awk '{if (NR==1){print}else{key=$1" "$2" "$3" "$5;if ($6>max[key]){max[key]=$6;line[key]=$0}}}END{for (k in line){print line[k]}}' he_eqtl.txt > he_eqtl.lead.txt
awk 'NR==1||($8>0&&$13<0)||($8<0&&$13>0)' he_eqtl.lead.txt > he_eqtl.lead.oppo.txt
