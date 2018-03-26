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
    blastadmin.py ln ${strain} ${inFilepath}
    blastadmin.py createdb blastn ${strain}
done

