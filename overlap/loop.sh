#!/bin/bash

waitfor() {
    filepath=${1}
    while true;
    do
        if [ -e ${filepath} ]
        then
            break
        fi
        sleep 1m
    done
}

myrun() {
    target=${1}
    statusFilename=.STATUS
    ./pipeline.sh ${target} ${statusFilename}
    
    echo "WAIT: ${statusFilename}"
    waitfor ${statusFilename}
    status=`cat ${statusFilename}`
    if [ ${status} -eq 0 ]; then
        echo "DONE: ${step}"
    else
        echo "ERROR: in ${step}"
    fi
    rm ${statusFilename}
}

while read line
do
    target=`echo ${line} | cut -d ',' -f1`
    if [ ${target:0:1} != ? ]; then
        echo ${target}
        myrun ${target}
    fi    
done < ../record.txt
