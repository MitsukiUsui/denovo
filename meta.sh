#!/bin/bash

MYVERSION="ver1.1"
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
#        echo "WAIT: for ${filepath}"
        sleep 1m
    done
}

echo "START: meta-pipeline for ${target} @${TIMESTAMP}"

#--------------------------------------------------------------------------------
# STEP0. download data
#--------------------------------------------------------------------------------
step=download
echo "START: ${step}"
cd ./${step}
./pipeline.sh ${target}
cd ../
echo "DONE: ${step}"

#--------------------------------------------------------------------------------
# STEP1. ortholog clustering
#--------------------------------------------------------------------------------
step=ortho
echo "START: ${step}"
cd ./${step}
./pipeline.sh ${target} ${statusFilename}

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

#--------------------------------------------------------------------------------
# STEP2. BLASTN
#--------------------------------------------------------------------------------
step=blastn
echo "START: ${step}"
cd ./${step}
./pipeline.sh ${target} ${statusFilename}

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

#--------------------------------------------------------------------------------
# STEP3. calculate synteny relationship
#--------------------------------------------------------------------------------
step=synteny
echo "START: ${step}"
cd ./${step}
./synteny.py ${target}
cd ../
echo "DONE: ${step}"

#--------------------------------------------------------------------------------
# STEP4. calculate DNA/Protein alignment score for overlap 
#--------------------------------------------------------------------------------
step=overlap
echo "START: ${step}"
cd ./${step}
./pipeline.sh ${target} ${statusFilename}

waitfor ${statusFilename}
status=`cat ${statusFilename}`
if [ ${status} -eq 0 ]; then
    echo "DONE ${step}"
    rm ${statusFilename}
    cd ../
else
    echo "ERROR in ${step}"
    exit ${status}
fi

echo "DONE: all steps successfully"
echo "${target},${MYVERSION},${TIMESTAMP},${DATE}" >> record.txt
