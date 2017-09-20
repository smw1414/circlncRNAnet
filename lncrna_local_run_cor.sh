#!/bin/bash
genes=$1
./correlation.r -n output/norm_readstable.txt -f count_matrix_conditon.txt -q $genes -a gencodev25 -s symbol -m lnc
./coexp_pairs_annotation.r -i output/linc_coexp_pairs.txt -o output/linc_coexp_pairs_ann.txt

