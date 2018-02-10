#!/bin/bash

IFS=$'\n'

target=${1}
outDirec=./db/${target}
mkdir -p ${outDirec}

echo "START: create databases"
catalogFilepath=../data/${target}/catalog.tsv
cut -f1 ${catalogFilepath} | tail -n +2 | while read strain
do
	inFilepath=/data/mitsuki/data/denovo/${target}/dnaseq/${strain}.dnaseq
	dbName=${outDirec}/${strain}
	makeblastdb -in ${inFilepath} -dbtype nucl -out ${dbName} -logfile /dev/null
	echo -e "\tDONE: output ${dbName}"
done

