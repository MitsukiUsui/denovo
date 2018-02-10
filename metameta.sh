#/bin/bash

while read target
do
    ./meta.sh ${target}
done < target.lst
