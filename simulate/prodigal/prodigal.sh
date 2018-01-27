#!/bin/bash

seqFilepath=${1}
gffFilepath=${2}
fnaFilepath=${3}
faaFilepath=${4}

time prodigal -i ${seqFilepath} \
              -o ${gffFilepath} -f gff \
              -d ${fnaFilepath} \
              -a ${faaFilepath}
