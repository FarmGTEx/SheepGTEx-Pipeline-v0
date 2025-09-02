# filter and concat vcf 
mkdir -p vcf
for chr in chr{1..26}
do
    jsub -q normal -n 4 -R "span[hosts=1]" -J vcf_${chr} \
        -e vcf/${chr}.%J.log -o vcf/${chr}.%J.log "bash vcf_filter.sh $chr 4"
done
jsub -q normal -n 4 -R "span[hosts=1]" -J vcf_concat \
    -e vcf/concat.%J.log -o vcf/concat.%J.log "bash vcf_concat.sh 4"

