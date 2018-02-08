#!/usr/bin/env python3

import pandas as pd
import os
import sys

sys.path.append("../helper")
from gff import read_gff, write_gff

def main(inFilepath, outFilepath):
    """
    add orf_id column to .gff created by prodigal
    """

    gff_df=read_gff(inFilepath, additional_lst=["ID"])
    attribute_lst=[]
    for _, row in gff_df.iterrows():
        orfId = "{}_{}".format(row["seqname"], row["ID"].split("_")[-1])
        att = "{};orf_id={}".format(row["attribute"], orfId)
        attribute_lst.append(att)
    gff_df["attribute"]=attribute_lst
    write_gff(gff_df, outFilepath)
    print("DONE: output {}".format(outFilepath))
    
if __name__=="__main__":
    gffFilepath=sys.argv[1]
    main(gffFilepath, gffFilepath)
