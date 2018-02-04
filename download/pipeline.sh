#!/bin/bash

target=${1}
./strain.sh ${target}
./arg/mywget.py ${target} > ./arg/mywget.lst
./mywget.sh ./arg/mywget.lst
