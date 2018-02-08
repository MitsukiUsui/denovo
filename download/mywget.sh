argFilepath=${1}
echo "START: wget & gunzip from refseq"
while read line
do
    ftpFilepath=`echo ${line} | cut -d ',' -f1`
    outFilepath=`echo ${line} | cut -d ',' -f2`
    wget --quiet -O - ${ftpFilepath} | gunzip -c > ${outFilepath}
    echo -e "\tDONE: ${outFilepath}"
done<${argFilepath}
