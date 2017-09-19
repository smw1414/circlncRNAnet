#!/bin/bash
rawread_count=$1
contab=$2
circ_count=$3
rm -rf output
mkdir -p output

./deseq2.r -r $rawread_count -f $contab -i $circ_count  -m circ -a gencodev25 -c 4
cp $contab count_matrix_conditon.txt
