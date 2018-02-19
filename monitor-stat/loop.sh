#!/bin/bash

while read line
do
    target=`echo ${line} | cut -d ',' -f1`
    if [ ${target:0:1} != ? ]; then
        echo ${target}
#        ./blastn_stat.py ${target}
        ./overlap_stat.py ${target}
    fi    
done < ../record.txt
