#!/bin/bash
#$ -S /bin/bash
#$ -N score
#$ -q all.q
#$ -cwd
#$ -v PATH
#$ -o ./log/$JOB_ID.out
#$ -e ./log/$JOB_ID.err
#$ -l mem_free=5G

cd ../
time ./find_overlap.py ${1} ${2}
time ./calc_identity.py ${1} ${2}
