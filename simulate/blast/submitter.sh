#!/bin/bash

train_array=("train2000" "train10-4" "train10-5" "train10-6")
dist_array=("0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9")
for trainName in ${train_array[@]}
do
    for dist in ${dist_array[@]}
    do
        seqName="sim+1_${dist}"
#        echo ${seqName},${trainName}
        ./create_database.sh ${seqName} ${trainName}
        ./blastn.sh ${seqName} ${trainName}
    done
done
