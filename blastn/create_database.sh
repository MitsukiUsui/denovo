#!/bin/bash

IFS=$'\n'

target=${1}
strainFilepath=../data/${target}/strain.lst
outDirec=./db/${target}
mkdir -p ${outDirec}

for strain in `cat ${strainFilepath}`
do
	inFilepath=/data/mitsuki/data/denovo/${target}/dnaseq/${strain}.dnaseq
	dbName=${outDirec}/${strain}
	makeblastdb -in ${inFilepath} -dbtype nucl -out ${dbName} -logfile /dev/null
	echo "CREATED ${dbName}"
done
