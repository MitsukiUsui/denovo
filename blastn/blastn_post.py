#!/usr/bin/env python3

import pandas as pd
import sys

def main(target, subStrain, lookupFilepath):
    lookup_df=pd.read_csv(lookupFilepath)
    
    hitFilepath="./result/{}/{}.tab".format(target, subStrain)
    columns_lst=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart","send", "evalue", "bitscore"]
    hit_df=pd.read_csv(hitFilepath, delimiter='\t', header=None, names=columns_lst)
    hit_df=hit_df.drop_duplicates()

    #add qfamily column
    tmp_df = lookup_df[~lookup_df[subStrain].isnull()][["family", subStrain]]
    tmp_df = tmp_df.rename(columns={"family": "qfamily", subStrain:"qseqid"})
    out_df = pd.merge(hit_df, tmp_df, on="qseqid", how="left")
    assert out_df["qfamily"].isnull().sum()==0

    #add sstrain, qstrain, hit_strand columns
    dct_lst=[]
    for _,row in out_df.iterrows():
        dct={}
        dct["sstrain"]=subStrain
        dct["qstrain"]=row["qseqid"].split('-')[0]
        if row["sstart"]<row["send"]:
            dct["hit_strand"]=1
        else:
            dct["hit_strand"]=-1
        dct_lst.append(dct)
    assert out_df.shape[0]==len(dct_lst)
    tmp_df=pd.DataFrame(dct_lst, index=out_df.index)
    out_df=pd.concat([out_df, tmp_df], axis=1)
    out_df.index.name="hit_id"

    outFilepath=hitFilepath.replace(".tab", ".csv")
    out_df.to_csv(outFilepath)
    print("DONE: output to {}".format(outFilepath))

if __name__=="__main__":
#    target="sulfolobus_islandicus"
    target=sys.argv[1]
    subStrain=sys.argv[2]
    lookupFilepath="../data/{}/query_lookup.csv".format(target)
    main(target, subStrain, lookupFilepath)
