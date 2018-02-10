#!/bin/bash

target=${1}
statusFilename=${2:-.STATUS}

./create_database.sh ${target}
./query_lookup.py ${target}
./create_query.py ${target}

cmd=blastn.sh
argCmd=./arg/${cmd%.*}.py
argFilepath=${argCmd%.*}.lst
eval ${argCmd} ${target} > ${argFilepath}
numJobs=`grep -c '' ${argFilepath}`
jobId_r=`qsub -terse -t 1-${numJobs} ${cmd} ${argFilepath}`
jobId=`echo ${jobId_r} | cut -d '.' -f1`
echo "submitted ${numJobs} jobs with job_id=${jobId}"

cmd=checker.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target} ${statusFilename}`
echo "submitted checker with job_id=${jobId}, dependency=${prevJobId}"

