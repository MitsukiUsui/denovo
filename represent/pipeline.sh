#!/bin/bash

target=${1}
statusFilename=${2:-.STATUS}

./myphylophlan.sh ${target}
./represent.py ${target} 0.5
echo 0 > ${statusFilename}
