#!/bin/bash
mkdir -p log splits
jsub -R "span[hosts=1]" -q normal -n 4 -J v4.03.store \
	-e log/03.store.%J.log -o log/03.store.%J.log \
	"bash 03.store.sh"
