#!/bin/bash

target=bacillus
inputDirec=/data/mitsuki/data/mbgd/dnaseq
seqDirec=/data/mitsuki/data/ortho/${target}/dnaseq
prodigalDirec=/data/mitsuki/data/ortho/${target}/prodigal

mkdir -p ${seqDirec}
mkdir -p ${prodigalDirec}

strainFilepath=/home/mitsuki/altorf/mbgd/data/${target}/strain.lst
while read strain
do
    inputFilepath=${inputDirec}/${strain}.dnaseq
    seqFilepath=${seqDirec}/${strain}.dnaseq
#    ln -s ${inputFilepath} ${seqFilepath}
    
    gffFilepath=${prodigalDirec}/${strain}.gff
    fnaFilepath=${prodigalDirec}/${strain}.fna
    faaFilepath=${prodigalDirec}/${strain}.faa
    supFilepath=${prodigalDirec}/${strain}.sup # supplementary information 
    qsub prodigal.sh ${seqFilepath} ${gffFilepath} ${fnaFilepath} ${faaFilepath} ${supFilepath}
done < ${strainFilepath}

