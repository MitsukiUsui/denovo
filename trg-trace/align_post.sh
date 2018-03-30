#!/bin/bash
#$ -S /bin/bash
#$ -N post
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/post_$JOB_ID.out
#$ -e ./log/post_$JOB_ID.err
#$ -l mem_free=5G

set -eu

target=${1}
./align_post.py ${target}

