#!/bin/bash
set -u

MYVERSION="ver2.0"
TIMESTAMP=`date +%s`
DATE=`date -d @${TIMESTAMP}`

target=${1}
statusFilename=.${TIMESTAMP}

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
    # require pipeline.sh & corresponting checker.sh mechanism
    
    step=${1}
   
    echo ""
    echo "START: ${step}"
    cd ./${step}
    ./pipeline.sh ${target} ${statusFilename}

    echo "WAIT: ${statusFilename}"
    waitfor ${statusFilename}
    status=`cat ${statusFilename}`
    if [ ${status} -eq 0 ]; then
        echo "DONE: ${step}"
        rm ${statusFilename}
        cd ../
    else
        echo "ERROR: in ${step}"
        exit ${status}
    fi
}

echo "START: meta-pipeline for ${target} @${TIMESTAMP}"

myrun download
myrun phylogeny
myrun ortho
myrun trg 
myrun blastn

echo "DONE: all steps successfully"
echo "${target},${MYVERSION},${TIMESTAMP},${DATE}" >> record.txt
