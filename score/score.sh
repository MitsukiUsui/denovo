#!/bin/bash
#$ -S /bin/bash
#$ -N score
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/score_$JOB_ID.out
#$ -e ./log/score_$JOB_ID.err
#$ -l mem_free=5G

target=${1}
./orf2score.py ${target}
./family2score.py ${target}
