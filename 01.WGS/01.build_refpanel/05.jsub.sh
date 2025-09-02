#!/bin/bash
jsub -q fat -M 70000000 -n 2 -R "span[hosts=1]" -J sheep.refpanel.snpEff \
        -e log/snpEff.%J.log -o log/snpEff.%J.log \
        "bash 05.snpEff.sh"

