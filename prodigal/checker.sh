#!/bin/bash
#$ -S /bin/bash
#$ -N cheker
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/checker_$JOB_ID.out
#$ -e ./log/checker_$JOB_ID.err
#$ -l mem_free=5G

target=${1}
statusFilename=${2}
./checker.py ${target}
echo $?> ${statusFilename}

