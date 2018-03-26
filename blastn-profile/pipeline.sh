#!/bin/bash

set -u

target=${1}
statusFilename=${2:-.STATUS}

cmd=region.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
jobid_r=`qsub -terse -t 1-${numJobs} ${cmd} ${argFilepath}`
jobid=`echo ${jobid_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobid}"

cmd=checker.sh
numJobs=1
prvJobid=${jobid}
jobid=`qsub -terse -hold_jid ${prvJobid} ${cmd} ${target} ${statusFilename}`
echo "submitted checker with job_id=${jobid}, dependency=${prvJobid}"

