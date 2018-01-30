#!/bin/bash

target=alteromonas
annotType=refseq

#--------------------------------------------------------------------------------
# STEP1. ortholog clustering
#--------------------------------------------------------------------------------
step=ortho
cd ${step}
./pipeline.sh ${target} ${annotType}
statusFilepath=.STATUS
rm ${statusFilepath}
while true
do
    if [ -e ${statusFilepath} ];
    then
        break
    fi
    echo "WAIT: for ${statusFilepath}"
	sleep 5s
done
status=`cat ${statusFilepath}`
if [ ${status} -eq 0 ]; then
    echo "DONE ${step}"
    cd ../
else
    echo "ERROR in ${step}"
    exit ${status}
fi
