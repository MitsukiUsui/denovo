#!/bin/bash

target=${1}
statusFilename=${2:-.STATUS}

#--------------------------------------------------------------------------------
# download from RefSeq-ftp and store
#--------------------------------------------------------------------------------
cmd=mywget.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
jobId_r=`qsub -terse -t 1-${numJobs} -tc 10 ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}"

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
prevJobId=${jobId}
jobId_r=`qsub -terse -hold_jid ${prevJobId} -t 1-${numJobs} ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}, dependency=${prevJobId}"

#--------------------------------------------------------------------------------
# checker
#--------------------------------------------------------------------------------
cmd=checker.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target} ${statusFilename}`
echo "submitted checker with job_id=${jobId}, dependency=${prevJobId}"
