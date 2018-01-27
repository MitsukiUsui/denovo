#!/bin/bash

target=${1}
inFilepath=./${target}/catalog.tsv
outFilepath=./${target}/strain.lst
awk -F "\t" '{if (NR > 1) print $1}' ${inFilepath} > ${outFilepath}
echo "DONE: output ${outFilepath}"
