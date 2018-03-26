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

argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
target=`echo ${line} | cut -d ',' -f1`
strain=`echo ${line} | cut -d ',' -f2`


dir=.
queryFilepath=${dir}/query/${target}/${strain}.fna
outFilepath=${dir}/result/${target}/${strain}.tab

mkdir -p `dirname ${outFilepath}`
blastadmin.py search blastn ${queryFilepath} ${strain} ${outFilepath}
./blastn_post.py ${target} ${strain}
