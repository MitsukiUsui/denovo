#!/bin/bash
#$ -S /bin/bash
#$ -N trg
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/trg_$JOB_ID.out
#$ -e ./log/trg_$JOB_ID.err
#$ -l mem_free=100G

set -u

target=${1}
./split.py ${target}
./split.sh ${target}
./lca.py ${target}
./query_lookup.py ${target}

