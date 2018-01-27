#!/bin/bash
#$ -S /bin/bash
#$ -N refseq
#$ -q all.q
#$ -cwd
#$ -v PATH
#$ -o ./log/$JOB_ID.out
#$ -e ./log/$JOB_ID.err
#$ -l mem_free=5G

cd ../
time ./format_refseq.py ${1} ${2} ${3}
