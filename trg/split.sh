#!/bin/bash

set -u
target=${1}

baseDirec=/data/mitsuki/out/altorf/denovo/trg/${target}
inFilepath=${baseDirec}/result.csv
outDirec=${baseDirec}/family

awk -F "," -v d="${outDirec}" 'NR>1{fn=sprintf("%s/%s.csv",d,$1); print $0 >> fn}' ${inFilepath}
