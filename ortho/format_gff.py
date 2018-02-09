#!/usr/bin/env python3

"""
add family column to .gff.
can be work when one orf_id belongs to multiple family, although this is not the case for sonicparanoid
"""

import pandas as pd
from io import StringIO
from collections import defaultdict
import os
import sys

sys.path.append("../helper")
from gff import read_gff, write_gff
from myio import get_strain_lst

def main(strain, clusterFilepath, inFilepath, outFilepath):
    cluster_df=pd.read_csv(clusterFilepath, delimiter="\t")
    filter_df=cluster_df[~cluster_df[strain].isnull()]
    orf2family=defaultdict(list) # key: orf_id , val: list of families which belongs to
    for family, orfIds in zip(filter_df["family"], filter_df[strain]):
        for orfId in orfIds.split(","):
            orf2family[orfId].append(family)
            
    # add family information to attribute
    gff_df=read_gff(inFilepath, additional_lst=["orf_id"])
    attribute_lst=[]
    assignCount = 0
    for _, row in gff_df.iterrows(): #!!! assuming every row is CDS
        if row["orf_id"] in orf2family.keys():
            assignCount += 1
            att = "{};family={}".format(row["attribute"], ",".join(orf2family[row["orf_id"]]) )
            attribute_lst.append(att)
        else: # when this CDS is pseudogene, no sequence information in FASTA so that no family information is given
#            print("DEBUG: {} has no family information".format(row["orf_id"]))
            attribute_lst.append(row["attribute"])
    assert assignCount == len(orf2family) #every orf_id should appear exactly once in gff, as orf_id is uniquely defined
    gff_df["attribute"]=attribute_lst
    write_gff(gff_df, outFilepath)
    print("DONE: output {}".format(outFilepath))
    print("\tassigned family to {}/{} CDS".format(assignCount, gff_df.shape[0]))
    
if __name__=="__main__":
    target=sys.argv[1]
    annType=sys.argv[2]
    clusterFilepath="../data/{}/cluster.tsv".format(target)    
    strain_lst = get_strain_lst(target)
    for strain in strain_lst:
        inFilepath="/data/mitsuki/data/denovo/{}/annotation/{}/gff/{}.gff".format(target, annType, strain)
        outFilepath = inFilepath 
        main(strain, clusterFilepath, inFilepath, outFilepath)
