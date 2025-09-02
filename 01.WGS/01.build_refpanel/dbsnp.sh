#!/bin/bash
# dbsnps downloaded from https://ftp.ensembl.org/pub/release-113/variation/gvf/ovis_aries_rambouillet
zcat ovis_aries_rambouillet.gvf.gz | awk '$3=="SNV"{print $1"\t"$4"\t"$9}' | sed -e 's/Variant_seq=/ /g' -e 's/:rs/ rs/g' -e 's/Reference_seq=/ /g' | awk '{print "chr"$1"\t"$2"\t"$1"_"$2"\t"$5"\t"$4"\t"$6}' | sed 's/;/ /g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$NF}' | sed '1iCHROM\tPOS\tvariant_id\trsid\tALT\tREF' | pigz -c > snpid2rsid.txt.gz
