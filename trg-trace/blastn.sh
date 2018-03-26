#!/bin/bash
#$ -S /bin/bash
#$ -N blastn
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -v BLASTADMIN_DATA
#$ -o ./log/blastn_$JOB_ID_$TASK_ID.out
#$ -e ./log/blastn_$JOB_ID_$TASK_ID.err
#$ -l mem_free=5G

set -ue

argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
id=`echo ${line} | cut -d ',' -f1`
filepath=`echo ${line} | cut -d ',' -f2`
query=`echo ${line} | cut -d ',' -f3`
result=`echo ${line} | cut -d ',' -f4`

yes no|blastadmin.py ln ${id} ${filepath}
blastadmin.py search blastn ${id} ${query} ${result}
