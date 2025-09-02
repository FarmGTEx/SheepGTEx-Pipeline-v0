#!/bin/bash
for j in `cut -f1 tissue.list`
do
	jsub -R "span[hosts=1]" -q normal -n 1 -J v4.02.eval.single.${j} \
		-e log/02.eval.single.${j}.%J.log -o log/02.eval.single.${j}.%J.log \
		"bash 02.eval.single.sh $j"
done
