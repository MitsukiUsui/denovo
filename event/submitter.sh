#!/bin/bash

while read target
do
    ./post_event.py ${target} &
done < target.lst
