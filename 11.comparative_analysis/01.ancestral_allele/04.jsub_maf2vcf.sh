for chrom in `cat chr.list`
do 
	for species in {'Cattle','Pig','Human'}
	do
		bash 04.maf2vcf.sh ${chrom} ${species}
	done
done
