#!/usr/bin/env python3

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
          
    # add orf_id & family information to attribute
    gff_df=read_gff(inFilepath, additional_lst=["ID"])
    attribute_lst=[]
    for _, row in gff_df.iterrows():
        orfId = "{}_{}".format(row["seqname"], row["ID"].split("_")[-1])
        att = "{};orf_id={};family={}".format(row["attribute"], orfId, orf2family[orfId])
        attribute_lst.append(att)
    gff_df["attribute"]=attribute_lst
    write_gff(outFilepath, gff_df)
    print("DONE: output {}".format(outFilepath))
    
if __name__=="__main__":
    target="bacillus"
    clusterFilepath="../data/{}/cluster.tsv".format(target)
    outDirec="/data/mitsuki/data/ortho/{}/gff".format(target)
    os.makedirs(outDirec, exist_ok=True)
    
    strainFilepath="../data/{}/strain.lst".format(target)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    for strain in strain_lst:
        inFilepath="/data/mitsuki/data/ortho/{}/prodigal/{}.gff".format(target, strain)
        outFilepath="{}/{}.gff".format(outDirec, strain)
        main(strain, clusterFilepath, inFilepath, outFilepath)
