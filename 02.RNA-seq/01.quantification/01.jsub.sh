#!/bin/bash
mkdir -p Logs
#sed 's/:/\t/g' errors/error.log | cut -f1 | xargs rm
job=$1
jsub -q fat -M 32000000 -n 1 -R "span[hosts=1]" -J sheep.MP.${job} \
		-e Logs/sheep.${job}.%J.log -o Logs/sheep.${job}.%J.log \
	       	"bash 02.command.sh ${job}"
