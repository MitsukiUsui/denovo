#!/bin/bash

seqName=${1}
seqFilepath=${seqName}.fna
alignFilepath=${1}.align.fna
mafft --auto ${seqFilepath} > ${alignFilepath}


