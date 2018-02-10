#!/bin/bash

target=${1}
./myphylophlan.sh ${target}
./represent.py ${target} 0.05
