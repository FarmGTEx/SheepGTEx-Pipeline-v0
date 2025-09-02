#!/bin/bash
cat list | while read species name version
do
	jsub -q normal -n 50 -R "span[hosts=1]" -J 02.lastal_sheep_${species}_${name} -e ${species}_${name}_${version}/lastal.%J.log -o ${species}_${name}_${version}/lastal.%J.log "bash 02.lastal.sh $species $name $version"
done
