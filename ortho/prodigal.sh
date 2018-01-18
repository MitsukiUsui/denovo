#!/bin/bash
#$ -S /bin/bash
#$ -N prodigal
#$ -q all.q
#$ -cwd
#$ -v PATH
#$ -o ./log/$JOB_ID.out
#$ -e ./log/$JOB_ID.err
#$ -l mem_free=5G

seqFilepath=${1}
gffFilepath=${2}
fnaFilepath=${3}
faaFilepath=${4}
supFilepath=${5}

time prodigal -i ${seqFilepath} \
              -o ${gffFilepath} -f gff \
              -d ${fnaFilepath} \
              -a ${faaFilepath} \
              -s ${supFilepath}
