#!/bin/bash

target=${1}
statusFilename=.STATUS
baseDirec=/data/mitsuki/data/denovo/${target}

#--------------------------------------------------------------------------------
# annotation by prodigal
#--------------------------------------------------------------------------------
mkdir -p ${baseDirec}/annotation/prodigal/gff
mkdir -p ${baseDirec}/annotation/prodigal/fna
mkdir -p ${baseDirec}/annotation/prodigal/faa
mkdir -p ${baseDirec}/annotation/prodigal/sup

cmd=prodigal.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
jobId_r=`qsub -terse -t 1-${numJobs} ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}"

#--------------------------------------------------------------------------------
# create SQLite database from sup
#--------------------------------------------------------------------------------
cmd=create_db.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
prevJobId=${jobId}
jobId_r=`qsub -terse -t 1-${numJobs} -hold_jid ${prevJobId} ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}, dependency=${prevJobId}"

cmd=checker.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target}`${statusFilepath}
echo "submitted ${numJobs} jobs with job_id=${jobId}, dependency=${prevJobId}"
