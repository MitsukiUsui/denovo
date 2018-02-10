#!/usr/bin/env python3

import sys
import re
import math
import pandas as pd
import numpy as np

sys.path.append("../helper")
from myio import *

def main(target, orfFilepath, familyFilepath):
    orf2score_df = pd.read_csv(orfFilepath)
    cluster_df = get_cluster_df(target)
    strain_lst = get_strain_lst(target)
    
    print("START: load {}".format(orfFilepath))
    orf2score={}
    for orf, score in zip(orf2score_df["orf_id"], orf2score_df["Total"]):
        orf2score[orf]=score
        
    print("START: aggregate scores for {} families each".format(cluster_df.shape[0]))
    dct_lst=[]
    for _, row in cluster_df.iterrows():
        dct = {"family": row["family"]}

        sumScore=0
        count=0
        for strain in strain_lst:
            if isinstance(row[strain], str):
                for orf in row[strain].split(','):
                    if not(math.isnan(orf2score[orf])):
                        sumScore += orf2score[orf]
                        count += 1
        if count > 0:
            dct["ave_score"] = sumScore / count
        else:
            dct["ave_score"] = np.nan
        dct["score_count"] = count
        dct_lst.append(dct)
        
    family2score_df = pd.DataFrame(dct_lst)
    family2score_df = family2score_df[["family", "ave_score", "score_count"]]
    family2score_df.to_csv(familyFilepath, index=False)
    print("DONE: output {}".format(familyFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    orfFilepath="../data/{}/orf2score.csv".format(target)
    familyFilepath="../data/{}/family2score.csv".format(target)
    main(target, orfFilepath, familyFilepath)
    
