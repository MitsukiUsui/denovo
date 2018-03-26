#!/bin/bash
#$ -S /bin/bash
#$ -N wget
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -v BLASTADMIN_DATA
#$ -o ./log/wget_$JOB_ID_$TASK_ID.out
#$ -e ./log/wget_$JOB_ID_$TASK_ID.err
#$ -l mem_free=5G

set -ue

argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
id=`echo ${line} | cut -d ',' -f1`
ftp=`echo ${line} | cut -d ',' -f2`

yes no|blastadmin.py wget ${id} ${ftp}
