#!/usr/bin/env python3

import sys
import os
import pandas as pd
from Bio import SeqIO

sys.path.append("../helper")
from myio import *

def output_query(target, lookupFilepath, outDirec, ext):
    full = False
    lookup_df = pd.read_csv(lookupFilepath)
    strain_lst = get_strain_lst(target, full)

    orf2record = {}
    for family in lookup_df["family"]:
        seqFilepath = "/data/mitsuki/data/denovo/{0}/annotation/refseq/family/{2}/{1}.{2}".format(target, family, ext)
        for rec in SeqIO.parse(seqFilepath, "fasta"):
            rec.description = family
            orf2record[rec.id] = rec

    for strain in strain_lst:
        outFilepath="{}/{}.{}".format(outDirec, strain, ext)
        with open(outFilepath, "w") as f:
            for orfid in lookup_df[strain].dropna():
                SeqIO.write(orf2record[orfid], f, "fasta")
    print("DONE: output {} {} queries under {}".format(len(strain_lst), ext, outDirec))


def main(target, lookupFilepath, outDirec):
    os.makedirs(outDirec, exist_ok=True)
    output_query(target, lookupFilepath, outDirec, "fna")
    output_query(target, lookupFilepath, outDirec, "faa")

if __name__=="__main__":
    target=sys.argv[1]
    lookupFilepath = "./lookup/{}.csv".format(target)
    outDirec = "./query/{}".format(target)
    main(target, lookupFilepath, outDirec)
