#!/bin/bash
#$ -S /bin/bash
#$ -N blastn 
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

OUT=./log/blastn_${JOB_ID}_${SGE_TASK_ID}.out
ERR=./log/blastn_${JOB_ID}_${SGE_TASK_ID}.err
dir=.
dbName=${dir}/db/${target}/${strain}
queryFilepath=${dir}/query/${target}/${strain}.fna
outFilepath=${dir}/result/${target}/${strain}.tab
mkdir -p `dirname ${outFilepath}`

#time blastn -db ${dbName}\
#            -query ${queryFilepath}\
#            -out ${outFilepath}\
#            -word_size 6\
#            -evalue 1e-3\
#            -outfmt 6 1>>${OUT} 2>>${ERR}
./blastn_post.py ${target} ${strain} 1>>${OUT} 2>>${ERR}
