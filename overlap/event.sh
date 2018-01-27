#!/bin/bash
#$ -S /bin/bash
#$ -N event
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/event_$JOB_ID.out
#$ -e ./log/event_$JOB_ID.err
#$ -l mem_free=5G

./event.py ${1}
