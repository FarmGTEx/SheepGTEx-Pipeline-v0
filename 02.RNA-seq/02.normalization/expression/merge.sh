#!/bin/bash
num=$1

# calculate individual TPM from featureCounts output
echo "calculate TPM"
cat $num | while read path indtis sample
do
	if [ ! -f "${path}.TPM" ]; then
		echo -e "calculate TPM of ${sample}"
		sed '2d' $path | csvtk -t -H uniq -f 1,2,3,4,5,6 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/$6*1000}' > ${path}.RPK
		sumrpk=`awk '{sum+=$7}END{print sum}' ${path}.RPK`
		awk -v sumrpk=$sumrpk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/sumrpk*1000000}' ${path}.RPK | sed "1iGeneid\tChr\tStart\tEnd\tStrand\tLength\t${sample}" > ${path}.TPM
	fi
done

# merge counts
echo "merge counts"
head -1 ${num} | while read path indtis sample ; do csvtk -t uniq -f 1,2,3,4,5,6 $path | csvtk rename -t -f "02.STAR2pass/${sample}/${sample}_Aligned.sortedByCoord.out.bam" -n "$indtis" > ${num}.count.txt ; done
sed '1d' ${num} | while read path indtis sample
do
	csvtk -t uniq -f 1,2,3,4,5,6 $path | csvtk rename -t -f "02.STAR2pass/${sample}/${sample}_Aligned.sortedByCoord.out.bam" -n "$indtis" > ${num}.count.${indtis}.txt
	csvtk join -t -f "1,2,3,4,5,6;1,2,3,4,5,6" -H -O --na 0 ${num}.count.txt ${num}.count.${indtis}.txt > ${num}.count.tmp.txt
	mv ${num}.count.tmp.txt ${num}.count.txt
	rm ${num}.count.${indtis}.txt
done

# merge TPM
echo "merge TPM"
head -1 ${num} | while read path indtis sample ; do csvtk -t uniq -f 1,2,3,4,5,6 ${path}.TPM | csvtk rename -t -f "${sample}" -n "$indtis" > ${num}.tpm.txt ; done
sed '1d' ${num} | while read path indtis sample
do
	csvtk -t uniq -f 1,2,3,4,5,6 ${path}.TPM | csvtk rename -t -f "${sample}" -n "$indtis" > ${num}.tpm.${indtis}.txt
	csvtk join -t -f "1,2,3,4,5,6;1,2,3,4,5,6" -H -O --na 0 ${num}.tpm.txt ${num}.tpm.${indtis}.txt > ${num}.tpm.tmp.txt
	mv ${num}.tpm.tmp.txt ${num}.tpm.txt
	rm ${num}.tpm.${indtis}.txt
done
