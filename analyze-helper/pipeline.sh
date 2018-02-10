#!/bin/bash

target=${1}
statusFilename=${2:-.STATUS}

./igv_gff.py ${target}
./score.sh ${target}
echo 0 > ${statusFilename}
