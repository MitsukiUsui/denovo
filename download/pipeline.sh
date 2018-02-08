#!/bin/bash

target=${1}
#./strain.sh ${target}

#--------------------------------------------------------------------------------
# download from RefSeq-ftp and store
#--------------------------------------------------------------------------------
#./arg/mywget.py ${target} > ./arg/mywget.lst
#./mywget.sh ./arg/mywget.lst

#--------------------------------------------------------------------------------
# format refseq annotation information
#--------------------------------------------------------------------------------
baseDirec=/data/mitsuki/data/denovo/${target}
mkdir -p ${baseDirec}/dnaseq
mkdir -p ${baseDirec}/annotation/refseq/gff
mkdir -p ${baseDirec}/annotation/refseq/fna
mkdir -p ${baseDirec}/annotation/refseq/faa

cmd=refseq.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
jobId_r=`qsub -terse -t 1-${numJobs} ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}"
