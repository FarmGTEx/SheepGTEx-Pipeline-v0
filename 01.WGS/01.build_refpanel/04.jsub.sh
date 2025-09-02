#!/bin/bash
jsub -q normal -M 500000000 -n 64 -R "span[hosts=1]" -J sheep.concat \
	-e log/concat.log -o log/concat.log bash 04.concat.sh 64
