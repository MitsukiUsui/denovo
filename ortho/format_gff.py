#!/usr/bin/env python3

"""
add family column to .gff.
"""

import pandas as pd
from io import StringIO
import os
import sys

sys.path.append("../helper")
from gff import read_gff, write_gff

def main(strain, clusterFilepath, inFilepath, outFilepath):
    cluster_df=pd.read_csv(clusterFilepath, delimiter="\t")
    filter_df=cluster_df[~cluster_df[strain].isnull()]
    orf2family={}
    for family, orfIds in zip(filter_df["family"], filter_df[strain]):
        for orfId in orfIds.split(","):
            orf2family[orfId]=family
          
    # add family information to attribute
    gff_df=read_gff(inFilepath, additional_lst=["orf_id"])
    attribute_lst=[]
    assignCount = 0
    for _, row in gff_df.iterrows(): #!!! assuming every row is CDS. need to fix later
        if row["orf_id"] in orf2family.keys():
            assignCount += 1
            att = "{};family={}".format(row["attribute"], orf2family[row["orf_id"]])
            attribute_lst.append(att)
        else:
            attribute_lst.append(row["attribute"])
    gff_df["attribute"]=attribute_lst
    write_gff(outFilepath, gff_df)
    print("DONE: output {}".format(outFilepath))
    print("\tassigned family to {}/{} CDS".format(assignCount, gff_df.shape[0]))
    
if __name__=="__main__":
    target=sys.argv[1]
    annType="prodigal"
    clusterFilepath="../data/{}/cluster.tsv".format(target)    
    strainFilepath="../data/{}/strain.lst".format(target)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    for strain in strain_lst:
        inFilepath="/data/mitsuki/data/denovo/{}/annotation/{}/gff/{}.gff".format(target, annType, strain)
        outFilepath = inFilepath 
        main(strain, clusterFilepath, inFilepath, outFilepath)
