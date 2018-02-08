target=${1}
dataDirec=../data/${target}
outFilepath=${dataDirec}/cluster.phb

#--------------------------------------------------------------------------------
# force mode
#--------------------------------------------------------------------------------
FORCE_MODE=false
forceFilepath=${outFilepath}
if [ "$FORCE_MODE" = false ] && [ -e ${forceFilepath} ]; then
    echo "PASS: target file already exists"
    exit
fi

# create symbolic link in input directory
phyloDirec=/home/mitsuki/software/phylophlan
mkdir -p ${phyloDirec}/input/${target}
strainFilepath=${dataDirec}/strain.lst
while read strain
do
    from=/data/mitsuki/data/denovo/${target}/annotation/refseq/faa/${strain}.faa
    to=${phyloDirec}/input/${target}/${strain}.faa
    rm ${to}
    ln -s ${from} ${to}
done < ${strainFilepath}

# run
cd ${phyloDirec}
./phylophlan.py -u ${target} --nproc 10 1>/dev/null 2>/dev/null

# create link for newick
nwkFilepath=${phyloDirec}/output/${target}/${target}.tree.nwk
ln -s ${nwkFilepath} ${outFilepath}
echo "DONE: ${outFilepath}"
