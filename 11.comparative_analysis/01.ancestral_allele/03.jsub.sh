#!/bin/bash
cat list | while read species name version
do
	jsub -q normal -n 2 -R "span[hosts=1]" -J sheep_chain_${species}_${name} -e ${species}_${name}_${version}/chain.%J.log -o ${species}_${name}_${version}/chain.%J.log "bash 03.chain.sh $species $name $version"
done
