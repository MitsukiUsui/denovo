#!/usr/bin/env python3

import sys
import pandas as pd
from collections import defaultdict

def ddct2dct(ddct, strain_lst):
    dct={"mmseqs_id": ddct["mmseqs_id"]}
    lineage, size=0, 0
    for strain in strain_lst:
        if len(ddct[strain])>0:
            lineage+=1
            size+=len(ddct[strain])
            dct[strain]=",".join(ddct[strain])
    dct["lineage"]=lineage
    dct["size"]=size
    return dct
    
def init_ddct(mmseqsId):
    ddct = defaultdict(list)
    ddct["mmseqs_id"]=mmseqsId
    return ddct

def main(strainFilepath, mmseqsFilepath, outFilepath):
    print("START: parse {}".format(mmseqsFilepath))
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    mmseqs_df=pd.read_csv(mmseqsFilepath, sep="\t", header=None)
    
    #create ddct(finally dct) for each cluster. key:strain, val:list of orf_id
    #value is list because multiple protain from same strain can be assign to the cluster
    dct_lst=[]
    prv=mmseqs_df[0].iloc[0]
    ddct=init_ddct(prv) 
    for _, row in mmseqs_df.iterrows():
        if row[0]!=prv:
            dct_lst.append(ddct2dct(ddct, strain_lst))
            prv=row[0]
            ddct=init_ddct(prv)
        strain=row[1].split(":")[0]
        ddct[strain].append(row[1])
    dct_lst.append(ddct2dct(ddct, strain_lst))
    
    out_df=pd.DataFrame(dct_lst)
    out_df["family"]=["family{}".format(i) for i in out_df.index]
    out_df=out_df[["family", "mmseqs_id", "lineage", "size"]+strain_lst]
    out_df.index.name="ClusterID"
    out_df.to_csv(outFilepath, sep="\t")
    print("DONE: output to {}".format(outFilepath))

if __name__=="__main__":
    strainFilepath=sys.argv[1]
    mmseqsFilepath=sys.argv[2]
    outFilepath=sys.argv[3]
    main(strainFilepath, mmseqsFilepath, outFilepath)
