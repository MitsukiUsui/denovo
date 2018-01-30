#!/bin/bash

target=${1}
baseDirec=/data/mitsuki/data/denovo/${target}

#--------------------------------------------------------------------------------
# format refseq annotation information
#--------------------------------------------------------------------------------
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

#--------------------------------------------------------------------------------
# annotation by prodigal
#--------------------------------------------------------------------------------
#mkdir -p ${baseDirec}/annotation/prodigal/gff
#mkdir -p ${baseDirec}/annotation/prodigal/fna
#mkdir -p ${baseDirec}/annotation/prodigal/faa
#
#cmd=prodigal.sh
#argCmd=./arg/${cmd%.*}.py
#argFilepath=${argCmd%.*}.lst
#eval ${argCmd} ${target} > ${argFilepath}
#numJobs=`grep -c '' ${argFilepath}`
#prevJobId=${jobId}
#jobId_r=`qsub -terse -t 1-${numJobs} -hold_jid ${prevJobId} ${cmd} ${argFilepath}`
#jobId=`echo ${jobId_r} | cut -d '.' -f1`
#echo "submitted ${numJobs} jobs with job_id=${jobId}, dependency=${prevJobId}"

#--------------------------------------------------------------------------------
# protein clustering 
#--------------------------------------------------------------------------------
cmd=sonic.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target}`
echo "submitted ${numJobs} jobs with job_id=${jobId}, dependency=${prevJobId}"

#--------------------------------------------------------------------------------
# check if pipeline above worked correctly
#--------------------------------------------------------------------------------
cmd=checker.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target}`
echo "submitted checker with job_id=${jobId}, dependency=${prevJobId}"

