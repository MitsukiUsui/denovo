#!/bin/bash

target=${1}
statusFilename=${2:-.STATUS}

./synteny.py ${target}
echo 0 > ${statusFilename}
