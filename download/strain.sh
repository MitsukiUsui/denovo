#!/bin/bash

target=${1}
inFilepath=../data/${target}/catalog.tsv
outFilepath=../data/${target}/strain.lst
awk -F "\t" '{if (NR > 1) print $1}' ${inFilepath} > ${outFilepath}
echo "DONE: output ${outFilepath}"
