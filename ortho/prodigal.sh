#!/bin/bash
#$ -S /bin/bash
#$ -N prodigal
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o /dev/null
#$ -e /dev/null
#$ -l mem_free=5G

argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
seqFilepath=`echo ${line} | cut -d ',' -f1`
gffFilepath=`echo ${line} | cut -d ',' -f2`
fnaFilepath=`echo ${line} | cut -d ',' -f3`
faaFilepath=`echo ${line} | cut -d ',' -f4`
supFilepath=`echo ${line} | cut -d ',' -f5`

OUT=./log/prodigal_${JOB_ID}_${SGE_TASK_ID}.out
ERR=./log/prodigal_${JOB_ID}_${SGE_TASK_ID}.err

prodigal -i ${seqFilepath} \
         -o ${gffFilepath} -f gff \
         -d ${fnaFilepath} \
         -a ${faaFilepath} \
         -s ${supFilepath} 1>>${OUT} 2>>${ERR}

./prodigal_post.py ${gffFilepath} 1>>${OUT} 2>>${ERR}

