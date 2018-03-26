#!/bin/bash
#$ -S /bin/bash
#$ -N mmseqs
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -v BLASTADMIN_DATA
#$ -o ./log/mmseqs_$JOB_ID.out
#$ -e ./log/mmseqs_$JOB_ID.err
#$ -l mem_free=100G
#$ -pe make 20

set -u

target=${1}

direc=/data/mitsuki/out/altorf/denovo/trg/${target}
queryFilepath=${direc}/mmseqs/query.faa
resultFilepath=${direc}/mmseqs/result.m8
outFilepath=${direc}/result.csv

./mmseqs_query.py ${target} ${queryFilepath}
blastadmin.py search mmseqs nr ${queryFilepath} ${resultFilepath}
./mmseqs_post.py ${target} ${queryFilepath} ${resultFilepath} ${outFilepath}
