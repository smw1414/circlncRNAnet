#!/bin/bash
gene=$1
./correlation.r -n output/norm_readstable.txt -f count_matrix_conditon.txt -q $gene -a gencodev25 -s symbol -m circ -c output/norm_readstable_circRNA.txt
