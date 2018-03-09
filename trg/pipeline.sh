#!/bin/bash

set -u

target=${1}
statusFilename=${2:-.STATUS}

./mmseqs_pre.py ${target}
./mmseqs.sh ${target}
./mmseqs_post.py ${target}

./split.py ${target}
./split.sh ${target}

./lca.py ${target}
./trg.py ${target}

#--------------------------------------------------------------------------------
# check if pipeline above worked correctly
#--------------------------------------------------------------------------------
#cmd=checker.sh
#numJobs=1
#prevJobId=${jobId}
#jobId=`qsub -terse -hold_jid ${prevJobId} ${cmd} ${target} ${statusFilename}`
#echo "submitted checker with job_id=${jobId}, dependency=${prevJobId}"

