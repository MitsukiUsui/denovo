#!/bin/bash
#$ -S /bin/bash
#$ -N refseq
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/refseq_$JOB_ID_$TASK_ID.out
#$ -e ./log/refseq_$JOB_ID_$TASK_ID.err
#$ -l mem_free=5G

argFilepath=${1}
lineNum=${SGE_TASK_ID:-1}
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
target=`echo ${line} | cut -d ',' -f1`
basename=`echo ${line} | cut -d ',' -f2`
genomeid=`echo ${line} | cut -d ',' -f3`

./refseq.py ${target} ${basename} ${genomeid}
