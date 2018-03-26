#!/bin/bash

target=${1}
statusFilename=${2:-.STATUS}

#--------------------------------------------------------------------------------
# protein clustering
#--------------------------------------------------------------------------------
cmd=sonic.sh
numJobs=1
jobId=`qsub -terse ${cmd} ${target}`
echo "submitted ${numJobs} jobs with job_id=${jobId}"

#--------------------------------------------------------------------------------
# check if pipeline above worked correctly
#--------------------------------------------------------------------------------
cmd=checker.sh
numJobs=1
prevJobId=${jobId}
jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target} ${statusFilename}`
echo "submitted checker with job_id=${jobId}, dependency=${prevJobId}"

