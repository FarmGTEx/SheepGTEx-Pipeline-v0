#!/bin/bash
mkdir -p log splits
jsub -R "span[hosts=1]" -q gpu -n 1 -gpgpu "1 mig=1" -J v4.01.train \
	-e log/01.train.%J.log -o log/01.train.%J.log \
	"bash 01.train.sh"
