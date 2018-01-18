#!/bin/bash
#$ -S /bin/bash
#$ -N blastn 
#$ -q all.q
#$ -cwd
#$ -v PATH
#$ -o ./log/$JOB_ID.out
#$ -e ./log/$JOB_ID.err
#$ -l mem_free=5G

cd ../
#time ./blastn.sh ${1} ${2}
time ./summarize_blastn.py ${1} ${2}
