#!/bin/bash
#$ -S /bin/bash
#$ -N cluster
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/cluster_$JOB_ID.out
#$ -e ./log/cluster_$JOB_ID.err
#$ -l mem_free=30G

target=${1}
annotationType=refseq
baseDirec=/data/mitsuki/data/denovo/${target}/annotation/${annotationType}
mmseqsDirec=${baseDirec}/mmseqs
mkdir -p ${mmseqsDirec}


# merge all the faa file into allFilepath
allFilepath=${mmseqsDirec}/${target}.faa
rm ${allFilepath}
strainFilepath=../data/${target}/strain.lst
while read strain
do
    faaFilepath=${baseDirec}/faa/${strain}.faa
    cat ${faaFilepath} >> ${allFilepath}
done < ${strainFilepath}
echo "DONE: merge faa to ${allFilepath}"

# run clustering with mmseqs
cd ${mmseqsDirec}
mmseqs createdb ${allFilepath} DB
mkdir tmp
mmseqs cluster DB clu tmp
mmseqs createseqfiledb DB clu clu_seq 
mmseqs result2flat DB DB clu_seq clu_seq.fasta
mmseqs createtsv DB DB clu clu.tsv
echo "DONE: calculate clu.tsv"
echo ""

# format clu.tsv to ortholog table
mmseqsFilepath=${mmseqsDirec}/clu.tsv
clusterFilepath=../data/${target}/cluster.tsv
cd -
./cluster_post.py ${strainFilepath} ${mmseqsFilepath} ${clusterFilepath}
echo "DONE: format ortholog table to ${clusterFilepath}"
./format_gff.py ${target}
echo "DONE: add family column to gff"

