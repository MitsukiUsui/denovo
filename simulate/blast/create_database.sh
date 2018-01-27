#!/bin/bash

seqName=${1}
trainName=${2}

inFilepath=./sbjct/${trainName}/${seqName}.fna
dbFilepath=./db/${trainName}/${seqName}

mkdir -p `dirname ${dbFilepath}`
makeblastdb -in ${inFilepath} -dbtype nucl -out ${dbFilepath}
echo "CREATED ${dbName}"
