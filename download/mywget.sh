argFilepath=${1}
while read line
do
    ftpFilepath=`echo ${line} | cut -d ',' -f1`
    outFilepath=`echo ${line} | cut -d ',' -f2`
    wget --quiet -O - ${ftpFilepath} | gunzip -c > ${outFilepath}
    echo "DONE: ${outFilepath}"
done<${argFilepath}
