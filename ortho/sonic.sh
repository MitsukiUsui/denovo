#!/bin/bash
#$ -S /bin/bash
#$ -N sonic
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/sonic_$JOB_ID.out
#$ -e ./log/sonic_$JOB_ID.err
#$ -l mem_free=30G
#$ -pe make 20

target=${1}
annotType=refseq
sonicDirec=/data/mitsuki/data/denovo/${target}/annotation/${annotType}/sonic

rm -rf ${sonicDirec}
./sonic_pre.py ${target} ${annotType}
sonicparanoid.py  -i ${sonicDirec}/in -o ${sonicDirec}/out -t 20
./sonic_post.py ${target} ${annotType}
./format_gff.py ${target} ${annotType}
./distribute_family.py ${target}
./igv_gff.py ${target}
