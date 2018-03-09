#!/bin/bash
#$ -S /bin/bash
#$ -N region
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/region_$JOB_ID_$TASK_ID.out
#$ -e ./log/region_$JOB_ID_$TASK_ID.err
#$ -l mem_free=5G

argFilepath=${1}
if [ -z ${SGE_TASK_ID+x} ] ; then lineNum=1; else lineNum=${SGE_TASK_ID}; fi;
line=`awk -v lineNum=$lineNum '{if (NR == lineNum) print $0}' ${argFilepath}`
target=`echo ${line} | cut -d ',' -f1`
strain=`echo ${line} | cut -d ',' -f2`

mkdir -p ./out/${target}/${strain}

./define_region.py ${target} ${strain}
./define_extract.py ${target} ${strain}
./extract_seq.py ${target} ${strain}
./assign_phase.py ${target} ${strain}
./profile.py ${target} ${strain}
