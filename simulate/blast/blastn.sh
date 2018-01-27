#!/bin/bash

seqName=${1}
trainName=${2}

queryFilepath="../out/query/${seqName}.fna"
dbFilepath="./db/${trainName}/${seqName}"
outFilepath="./result/${trainName}/${seqName}.tab"
mkdir -p `dirname ${outFilepath}`

blastn -db ${dbFilepath}\
       -query ${queryFilepath}\
       -out ${outFilepath}\
       -word_size 6\
       -evalue 1e-1\
       -outfmt 6
echo "DONE: ${outFilepath}"
