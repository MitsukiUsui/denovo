#!/bin/bash

name=${1}
seqFilepath=./data/dnaseq/${name}.dnaseq
gffFilepath=./data/gff/${name}.gff
fnaFilepath=./data/fna/${name}.fna
faaFilepath=./data/faa/${name}.faa

tmpSeqFilepath=./prodigal/`basename ${seqFilepath}`
tmpGffFilepath=./prodigal/`basename ${gffFilepath}`
tmpFnaFilepath=./prodigal/`basename ${fnaFilepath}`
tmpFaaFilepath=./prodigal/`basename ${faaFilepath}`

#merge with train file
trainFilepath=./prodigal/train.dnaseq
cp ${trainFilepath} ${tmpSeqFilepath}
cat ${seqFilepath} >> ${tmpSeqFilepath}

prodigal -i ${tmpSeqFilepath} \
         -o ${tmpGffFilepath} -f gff \
         -d ${tmpFnaFilepath} \
         -a ${tmpFaaFilepath}

#filter result
./prodigal/filter_gff.py ${tmpGffFilepath} ${gffFilepath}
echo "DONE: ${gffFilepath}"
cat ${tmpFnaFilepath}| seqkit grep -rvp train > ${fnaFilepath}
echo "DONE: ${fnaFilepath}"
cat ${tmpFaaFilepath}| seqkit grep -rvp train > ${faaFilepath}
echo "DONE: ${faaFilepath}"
