#!/bin/bash
rawread_count=$1
contab=$2
rm -rf output
mkdir -p output
./deseq2.r -r $rawread_count -f $contab -m lnc -a gencodev25 -c 4
./deseq2_annotation.r -i output/DEGs_lncRNA.txt -o output/DEGs_lncRNA_ann.txt
cp $contab count_matrix_conditon.txt




