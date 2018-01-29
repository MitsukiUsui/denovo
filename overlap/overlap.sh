#!/bin/bash
#$ -S /bin/bash
#$ -N overlap
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o /dev/null
#$ -e /dev/null
#$ -l mem_free=5G

argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
target=`echo ${line} | cut -d ',' -f1`
strain=`echo ${line} | cut -d ',' -f2`

OUT=./log/overlap_${JOB_ID}_${SGE_TASK_ID}.out
ERR=./log/overlap_${JOB_ID}_${SGE_TASK_ID}.err
#./find_overlap.py ${target} ${strain} 1>>${OUT} 2>>${ERR}
#./calc_identity.py ${target} ${strain} 1>>${OUT} 2>>${ERR}
./assign_phase.py ${target} ${strain} 1>>${OUT} 2>>${ERR}
