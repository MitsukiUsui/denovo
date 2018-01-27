#!/bin/bash

strain=${1}

fnaFilepath=./${strain}.fna
cat ${fnaFilepath} | seqkit grep -rp ${strain}:sim > tmp
mv tmp ${fnaFilepath}

faaFilepath=./${strain}.faa
cat ${faaFilepath} | seqkit grep -rp ${strain}:sim > tmp
mv tmp ${faaFilepath}
