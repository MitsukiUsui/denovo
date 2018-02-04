#!/bin/bash

target=${1}
annotType=refseq
statusFilepath=.STATUS

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

echo "START: meta-pipeline for ${target},${annotType}"

#step=download
#echo "START: ${step}"
#cd ./${step}
#./pipeline.sh ${target}
#../
#echo "DONE: ${step}"

#--------------------------------------------------------------------------------
# STEP1. ortholog clustering
#--------------------------------------------------------------------------------
step=ortho
echo "START: ${step}"
cd ./${step}
rm ${statusFilepath}
./pipeline.sh ${target} ${annotType}

waitfor ${statusFilepath}
status=`cat ${statusFilepath}`
if [ ${status} -eq 0 ]; then
    echo "DONE: ${step}"
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
rm ${statusFilepath}
./pipeline.sh ${target} ${annotType}

waitfor ${statusFilepath}
status=`cat ${statusFilepath}`
if [ ${status} -eq 0 ]; then
    echo "DONE: ${step}"
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
rm ${statusFilepath}
./pipeline.sh ${target}

waitfor ${statusFilepath}
status=`cat ${statusFilepath}`
if [ ${status} -eq 0 ]; then
    echo "DONE ${step}"
    cd ../
else
    echo "ERROR in ${step}"
    exit ${status}
fi

echo "DONE: all steps successfully"
