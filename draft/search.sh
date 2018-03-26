#!/bin/bash
#$ -S /bin/bash
#$ -N search
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -v BLASTADMIN_DATA
#$ -o ./log/search_$JOB_ID_$TASK_ID.out
#$ -e ./log/search_$JOB_ID_$TASK_ID.err
#$ -l mem_free=5G

set -ue

argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
id=`echo ${line} | cut -d ',' -f1`
query=`echo ${line} | cut -d ',' -f2`
result=`echo ${line} | cut -d ',' -f3`

blastadmin.py search blastn ${id} ${query} ${result}
