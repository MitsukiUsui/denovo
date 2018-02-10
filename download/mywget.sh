#!/bin/bash
#$ -S /bin/bash
#$ -N wget
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/wget_$JOB_ID.out
#$ -e ./log/wget_$JOB_ID.err
#$ -l mem_free=5G

argFilepath=${1}
lineNum=${SGE_TASK_ID:-1}
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
ftpFilepath=`echo ${line} | cut -d ',' -f1`
outFilepath=`echo ${line} | cut -d ',' -f2`

wget --timeout 120 --no-verbose -O - ${ftpFilepath} | gunzip -c > ${outFilepath}

