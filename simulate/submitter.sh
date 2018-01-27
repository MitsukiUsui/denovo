#!/bin/bash

while read line
do
    ./prodigal_wrapper.sh ${line}
done < target.lst
