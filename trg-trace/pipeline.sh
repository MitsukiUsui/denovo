#!/bin/bash

set -ue

target=${1}

./create_query.py ${target}

mkdir -p ./result/${target}

cmd=blastn.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
jobId_r=`qsub -terse -t 1-${numJobs} ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}"
