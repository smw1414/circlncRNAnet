#!/bin/bash
tcga=$1
rm -rf output
mkdir -p output
./deseq2.r -m $tcga 
./deseq2_annotation.r -i output/DEGs_lncRNA.txt -o output/DEGs_lncRNA_ann.txt


