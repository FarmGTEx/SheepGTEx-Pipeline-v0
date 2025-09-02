for i in `cat chr.list` ;do gunzip -c chr${i}_ancestral_state.tsv.gz >> merged_ancestral_state.tsv; done
awk '!/^Chrom/ {print $0 "\t" $1"_"$2}' merged_ancestral_state.tsv > sheep_ancestral_state.tsv
rm -f merged_ancestral_state.tsv
awk 'NR==FNR{a[$2]; next} {if ($11 in a) print $0}' chrAuto.filtered.bim sheep_ancestral_state.tsv > sheep_overlap_ancestral.tsv
awk '{if ($9 <0.2) {count_low++;} else if ($9 >0.8) {count_high++;} else if ($9 >=0.2 && $9 <=0.8) {count_mid++;}} END {print "Minor allele <0.2: " count_low; print "Major allele>0.8: " count_high; print "Ambiguous allele >=0.2 and <=0.8: " count_mid;}' sheep_overlap_ancestral.tsv>sheep_summery_overlap.tsv
