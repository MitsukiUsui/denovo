#!/usr/bin/env python3

import sys
import pandas as pd
from ete3 import Tree
import numpy as np
from event_stat import EventStat

def get_family2score(target):    
    fp = "../data/{}/family2score.csv".format(target)
    df = pd.read_csv(fp)
    
    family2score = {}
    for family, score in zip(df["family"], df["ave_score"]):
        family2score[family]=score
    return family2score

def main(target, inFilepath, outFilepath):
    event_df = pd.read_csv(inFilepath, dtype = {"phase": object})
    print("DONE: load {} event from {}".format(event_df.shape[0], inFilepath))
    
    eventStat = EventStat(target)
    family2score= get_family2score(target)
    
    dct_lst = []
    for _, row in event_df.iterrows():
        dct = {"event_id": row["event_id"]}
        dct["family_distance"] = eventStat.family_distance(row["left"], row["right"])
        dct["left_score"] = family2score[row["left"]]
        dct["right_score"] = family2score[row["right"]]
        dct_lst.append(dct)
    tmp_df=pd.DataFrame(dct_lst)
    event_df = pd.merge(event_df, tmp_df, on = "event_id")
    
    event_df.to_csv(outFilepath, index = False)
    print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    inFilepath = "../overlap/out/{}/event.csv".format(target)
    outFilepath = "./out/{}.csv".format(target)
    main(target, inFilepath, outFilepath)