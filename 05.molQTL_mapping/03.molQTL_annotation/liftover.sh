#!/bin/bash
chunk=$1

liftOver -minMatch=0.2 -minBlocks=0.5 ${chunk} download/hg38ToGCF_016772045.2.over.chain.gz ${chunk}.Mapped ${chunk}.unMapped
pigz ${chunk}.Mapped
pigz ${chunk}.unMapped
