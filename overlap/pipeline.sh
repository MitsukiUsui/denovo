#!/bin/bash

target=${1}
statusFilename=${2:-.STATUS}

cmd=overlap.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
jobId_r=`qsub -terse -t 1-${numJobs} ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}"

#--------------------------------------------------------------------------------
# summarize into event information
#--------------------------------------------------------------------------------
cmd=event.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target}`
echo "submitted ${numJobs} jobs with job_id=${jobId}, dependency=${prevJobId}"

cmd=checker.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target} ${statusFilename}`
echo "submitted checker with job_id=${jobId}, dependency=${prevJobId}"
