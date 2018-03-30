#!/usr/bin/env python3

import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sys.path.append("../helper")
from myio import get_strain_lst

def main(target, strain, outDirec):

    # filter trace result
    traceFilepath = "./out/{}/{}.trace".format(target, strain)
    trace_df = pd.read_csv(traceFilepath)
    filtered_df = trace_df[trace_df["traceable"]==1]

    # create id2rec
    id2rec = {}
    queryFilepath = "./query/{}/{}.fna".format(target, strain)
    seqFilepath = "/data/mitsuki/data/denovo/{}/dnaseq/{}.dnaseq".format(target, strain)
    for fp in (queryFilepath, seqFilepath):
        for rec in SeqIO.parse(fp, "fasta"):
            id2rec[rec.id] = rec

    for _, row in filtered_df.iterrows():
        #extract query
        query = id2rec[row["query"]]

        #extract sbjct
        loci = row["longest"].split(';')
        assert len(loci) == 4
        id_ = "{}:[{}, {}),{}".format(loci[0], loci[1], loci[2], loci[3])
        seq = id2rec[loci[0]].seq[int(loci[1]) : int(loci[2])]
        if loci[3] == "1":
            seq = seq.reverse_complement()
        sbjct = SeqRecord(seq, id = id_, description = "")

        #output query, sbjct
        outFilepath = "{}/{}.fna".format(outDirec, row["query"])
        with open(outFilepath, "w") as f:
            SeqIO.write(query, f, "fasta")
            SeqIO.write(sbjct, f, "fasta")
#        print("DONE: output {}".format(outFilepath))
    print("DONE: prepare for total {}/{} alignments under {}".format(filtered_df.shape[0], trace_df.shape[0], outDirec))

if __name__=="__main__":
    target = sys.argv[1]
    strain_lst = get_strain_lst(target)
    for strain in strain_lst:
        outDirec = "./align/{}/{}".format(target, strain)
        os.makedirs(outDirec, exist_ok=True)
        main(target, strain, outDirec)
