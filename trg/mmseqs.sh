#!/bin/bash
#$ -S /bin/bash
#$ -N mmseqs
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/mmseqs_$JOB_ID.out
#$ -e ./log/mmseqs_$JOB_ID.err
#$ -l mem_free=100G
#$ -pe make 20

set -u

target=${1}

./mmseqs_pre.py ${target}

mmseqsDirec=/data/mitsuki/out/altorf/denovo/trg/${target}/mmseqs
seqFilepath=${mmseqsDirec}/query.faa
queryDB=${mmseqsDirec}/queryDB
resultDB=${mmseqsDirec}/resultDB
resultTsv=${mmseqsDirec}/result.m8
tmpDirec=${mmseqsDirec}/tmp
targetDB=/data/mitsuki/data/refseq/nr/targetDB

mkdir -p ${tmpDirec}
mmseqs createdb ${seqFilepath} ${queryDB}
mmseqs search ${queryDB} ${targetDB} ${resultDB} ${tmpDirec} --threads 20
mmseqs convertalis ${queryDB} ${targetDB} ${resultDB} ${resultTsv} --threads 20

./mmseqs_post.py ${target}
