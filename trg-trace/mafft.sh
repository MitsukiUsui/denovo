#!/bin/bash
#$ -S /bin/bash
#$ -N mafft
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o /dev/null
#$ -e /dev/null
#$ -l mem_free=5G

set -eu

#--------------------------------------------------------------------------------
# array job expander
#--------------------------------------------------------------------------------
argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
fnaFilepath=`echo ${line} | cut -d ',' -f1`
alignFilepath=`echo ${line} | cut -d ',' -f2`

mafft --auto ${fnaFilepath} > ${alignFilepath}
