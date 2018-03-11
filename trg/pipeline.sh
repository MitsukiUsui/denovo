#!/bin/bash

set -u

target=${1}
statusFilename=${2:-.STATUS}

cmd=mmseqs.sh
numJobs=1
jobid=`qsub -terse ${cmd} ${target}`
echo "submitted ${numJobs} jobs with job_id=${jobid}"

cmd=trg.sh
numJobs=1
prvJobid=${jobid}
jobid=`qsub -terse -hold_jid ${prvJobid} ${cmd} ${target}`
echo "submitted checker with job_id=${jobid}, dependency=${prvJobid}"

cmd=checker.sh
numJobs=1
prvJobid=${jobid}
jobid=`qsub -terse -hold_jid ${prvJobid} ${cmd} ${target} ${statusFilename}`
echo "submitted checker with job_id=${jobid}, dependency=${prvJobid}"

