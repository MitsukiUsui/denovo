#!/bin/bash
#$ -S /bin/bash
#$ -N blastn 
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/$JOB_ID.out
#$ -e ./log/$JOB_ID.err
#$ -l mem_free=5G

target=${1}
strain=${2}

dir=/home/mitsuki/altorf/denovo/blastn
dbName=${dir}/db/${target}/${strain}
queryFilepath=${dir}/query/${target}/${strain}.query
outFilepath=${dir}/result/${target}/${strain}.tab

mkdir -p `dirname ${outFilepath}`
blastn -db ${dbName}\
	   -query ${queryFilepath}\
	   -out ${outFilepath}\
	   -word_size 6\
	   -evalue 1e-3\
	   -outfmt 6 

