target=${1}

dataDirec=/home/mitsuki/altorf/denovo/data/${target}
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

catalogFilepath=${dataDirec}/catalog.tsv
cut -f1 ${catalogFilepath} | tail -n +2 | while read strain
do
    from=/data/mitsuki/data/denovo/${target}/annotation/refseq/faa/${strain}.faa
    to=${phyloDirec}/input/${target}/${strain}.faa
    rm ${to}
    ln -s ${from} ${to}
done

# run
cd ${phyloDirec}
./phylophlan.py -u ${target} --nproc 10 #1>/dev/null 2>/dev/null

# create link for newick
nwkFilepath=${phyloDirec}/output/${target}/${target}.tree.nwk
ln -s ${nwkFilepath} ${outFilepath}
echo "DONE: ${outFilepath}"
