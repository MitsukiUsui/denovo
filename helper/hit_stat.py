#!/usr/bin/env python3

import pandas as pd
import sys

def hit_stat(target, strain, family):
    baseDirec = "/home/mitsuki/altorf/denovo"
    lookupFilepath="{}/data/{}/query_lookup.csv".format(baseDirec, target)
    lookup_df=pd.read_csv(lookupFilepath)
    lookup_df = lookup_df[lookup_df["family"] == family]
    if lookup_df.shape[0] == 0:
        print("ERROR: no query for {}".format(family), file = sys.stderr)
        return
    else:
        assert lookup_df.shape[0] == 1
        query = lookup_df[lookup_df["family"] == family][strain].iloc[0]
        if not(isinstance(query, str)):
            print("ERROR: {} has {}".format(strain, family), file = sys.stderr)
            return
            
    print("Query: {}".format(query))
    hit_df=pd.read_csv("{}/blastn/result/{}/{}.csv".format(baseDirec, target, strain))
    hit_df = hit_df[hit_df["qfamily"] == family]
    ovr_df = pd.read_csv("{}/overlap/out/{}/{}_ovr.csv".format(baseDirec, target, strain), dtype = {"relative": object})
    ovr_df["score_dna"] = 0
    ovr_df["score_pro"] = 0
    col_lst = ["region_id", "ostart", "oend", "olength", "sfamily", "score_dna", "score_pro", "relative"]
    join_df = pd.merge(hit_df[hit_df["qfamily"]==family], ovr_df[col_lst], on = "region_id", how = "inner")
    
    print("{} HIT found:".format(hit_df.shape[0]))
    for _, row in hit_df.iterrows():
        stt = min(row["sstart"], row["send"])
        end = max(row["sstart"], row["send"])
        print("\t{:7d}:[{:7d}, {:7d})".format(row["length"], stt, end))

    print("{} OVR fount:".format(join_df.shape[0]))
    for _, row in join_df.iterrows():
        print("\t{:<11}, {:7d}:[{:7d},{:7d}), score:({:04.2f}, {:04.2f}), {}".format(row["sfamily"],  row["olength"], row["ostart"], row["oend"], row["score_dna"], row["score_pro"], row["relative"]))
        
if __name__=="__main__":
    hit_stat(sys.argv[1], sys.argv[2], sys.argv[3])
    
