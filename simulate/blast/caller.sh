#!/bin/bash

./create_database.sh ${1} ${2}
./blastn.sh ${1} ${2}
