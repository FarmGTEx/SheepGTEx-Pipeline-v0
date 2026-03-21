#!/bin/bash
### Usage: Merge RNA-Seq fastq file with cat. If a sample contains both single-end and pair-end data, we will only remain pair-end data. If a file is no need to merge, the file will be linked. 
### 
# help message
if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    sed -rn 's/^### ?//;T;p' "$0"
    echo """Usage: bash $0 <samplefile> <outdir> <sample>
Options:
    <samplefile>  A tsv file of all sample information that contains at most three columns (no header): Sample name, ablolute path of the sample fastq file of read1, ablolute path of the sample fastq file of read2 (no information for single-end data). One per line.
    <outdir>      The directory to store output sample path and files.
    <sample>      Name of the sample that will be processed.
    -h            Show this message.
    """
    exit 1
fi

# variable definition
samplefile=$1
outdir=$2
sample=$3

function mergefq() {
    mkdir -p ${outdir}/${sample}
    linenum=`grep -w "^$sample" $samplefile | wc -l`
    if [ $linenum -gt 1 ];then \
        # multiple runs
        pelinenum=`grep -w "^$sample" $samplefile | awk 'NF==3' | wc -l`
        if [ $pelinenum -eq 1 ];then \
            # one pair-end run with single-end run(s). Only remain the pair-end run.
            fq1=`grep -w "^$sample" $samplefile | awk 'NF==3{print $2}'`
            fq2=`grep -w "^$sample" $samplefile | awk 'NF==3{print $3}'`
            cmd1="ln -s $fq1 ${outdir}/${sample}/${sample}_1.fq.gz"
            echo $cmd1 ; $cmd1
            cmd2="ln -s $fq2 ${outdir}/${sample}/${sample}_2.fq.gz"
            echo $cmd2 ; $cmd2
        elif [ $pelinenum -gt 1 ];then \
            # multiple pair-end runs, merge them and ignore single-end run(s).
            echo "$sample read1 need to be merged:"
            echo `grep -w "^$sample" $samplefile | awk 'NF==3{print $2}'`
            grep -w "^$sample" $samplefile | awk 'NF==3{print $2}' | xargs cat > ${outdir}/${sample}/${sample}_1.fq.gz
            echo "$sample read2 need to be merged:"
            echo `grep -w "^$sample" $samplefile | awk 'NF==3{print $3}'`
            grep -w "^$sample" $samplefile | awk 'NF==3{print $3}' | xargs cat > ${outdir}/${sample}/${sample}_2.fq.gz
        else
            # multiple single-end runs, merge.
            echo "$sample read need to be merged:"
            echo `grep -w "^$sample" $samplefile | awk '{print $2}'`
            grep -w "^$sample" $samplefile | awk '{print $2}' | xargs cat > ${outdir}/${sample}/${sample}.fq.gz
        fi
        md5sum ${outdir}/${sample}/${sample}*.gz > ${outdir}/${sample}/md5.txt
    elif [ $linenum -eq 1 ];then \
        # single run
        fq1=`grep -w "^$sample" $samplefile | awk '{print $2}'`
        fq2=`grep -w "^$sample" $samplefile | awk '{print $3}'`
        if [ -f "$fq2" ];then \
            # pair-end
            cmd1="ln -s $fq1 ${outdir}/${sample}/${sample}_1.fq.gz"
            echo $cmd1 ; $cmd1
            cmd2="ln -s $fq2 ${outdir}/${sample}/${sample}_2.fq.gz"
            echo $cmd2 ; $cmd2
        else
            # single-end
            cmd="ln -s $fq1 ${outdir}/${sample}/${sample}.fq.gz"
            echo $cmd ; $cmd
        fi
        md5sum ${outdir}/${sample}/${sample}*.gz > ${outdir}/${sample}/md5.txt
    else
        echo "FAILED. No file for $sample, please check."
        exit 1
    fi
}

# main processes
function main() {
    mergefq
}

main
