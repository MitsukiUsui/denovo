#!/usr/bin/env python3

import sys

sys.path.append("../helper")
from gff import read_gff, write_gff

def filter_gff(inFilepath, outFilepath):
    gff_df=read_gff(inFilepath)
    filtered_df=gff_df[gff_df["seqname"] != "train"]
    write_gff(outFilepath, filtered_df)

if __name__=="__main__":
	inFilepath=sys.argv[1]
	outFilepath=sys.argv[2]
	filter_gff(inFilepath, outFilepath)

