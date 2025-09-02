#!/bin/bash
jsub -q normal -n 24 -R "span[hosts=1]" -J 01.lastdb.sh -e log/01.lastdb.sh.%J.log -o log/01.lastdb.sh.%J.log "bash 01.lastdb.sh"
