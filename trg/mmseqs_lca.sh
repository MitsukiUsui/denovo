#!/bin/bash
#$ -S /bin/bash
#$ -N mmseqs
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/mmseqs.out
#$ -e ./log/mmseqs.err
#$ -l mem_free=100G
#$ -pe make 20

set -u

target=nocardiaceae
seqFilepath=/data/mitsuki/data/denovo/${target}/annotation/refseq/faa/${target}.faa
mmseqsDirec=./out/lca/${target}
queryDB=${mmseqsDirec}/queryDB
targetDB=/data/mitsuki/data/uniprot/uniref90/targetDB
taxDirec=/data/mitsuki/data/refseq/taxonomy
lcaDB=${mmseqsDirec}/lcaDB
lcaTsv=${mmseqsDirec}/lca.tsv
tmpDirec=${mmseqsDirec}/tmp

mkdir -p ${mmseqsDirec}
mkdir -p ${tmpDirec}

#mmseqs createdb ${seqFilepath} ${queryDB}
mmseqs taxonomy ${queryDB} ${targetDB} ${targetDB}.tsv ${taxDirec} ${lcaDB} ${tmpDirec} --threads 30
mmseqs createtsv ${queryDB} ${lcaDB} ${lcaTsv}

