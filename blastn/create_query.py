#!/usr/bin/env python3

import sys
import os
import pandas as pd
from Bio import SeqIO

def output_query(ext, strain_lst, lookup_df):
    annotationType = "prodigal"

    print("START: load {} {} seqs".format(len(strain_lst), ext))
    orf2record={}
    for strain in strain_lst:
        seqFilepath="/data/mitsuki/data/denovo/{0}/annotation/{1}/{2}/{3}.{2}".format(target, annotationType, ext, strain)
        for record in SeqIO.parse(seqFilepath, "fasta"):
            orfId = record.id
            orf2record[orfId]=record
    
    print("START: create queries")
    for strain in strain_lst:
        outFilepath="./query/{}/{}.{}".format(target, strain, ext)
        with open(outFilepath, "w") as f:
            for orfId in lookup_df[strain].dropna():
                SeqIO.write(orf2record[orfId], f, "fasta")
        print("\tDONE: output {}".format(outFilepath))
    

def main(target, lookupFilepath, strainFilepath):
    lookup_df=pd.read_csv(lookupFilepath)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    
    outDirec = "./query/{}".format(target)
    os.makedirs(outDirec, exist_ok=True)
    output_query("fna", strain_lst, lookup_df)
    output_query("faa", strain_lst, lookup_df)

if __name__=="__main__":
    target=sys.argv[1]
    lookupFilepath="../data/{}/query_lookup.csv".format(target)
    strainFilepath="../data/{}/strain.lst".format(target)
    main(target, lookupFilepath, strainFilepath)
