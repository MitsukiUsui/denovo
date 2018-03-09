#!/usr/bin/env python3

import pandas as pd
import sys

sys.path.append("../helper")
from myio import *

def main(target, trgFilepath, outFilepath):
    strain_lst = get_strain_lst(target)
    cluster_df = get_cluster_df(target)
    trg_df = pd.read_csv(trgFilepath)
    cluster_df = pd.merge(cluster_df, trg_df, on="family", how="left")
    
    # create query for every member of trg (highly redundant though)
    dct_lst=[]
    for _, row in cluster_df.iterrows():
        if row["genus"] > 0 and row["lineage"] > 1:
            orfids_lst = list(row[strain_lst].dropna().values)
            for orfids in orfids_lst:
                for orfid in orfids.split(','):
                    dct={ "family" : row["family"]}
                    for strain in row[strain_lst].index[row[strain_lst].isnull()]:
                        dct[strain] = orfid
                    dct_lst.append(dct)
    
    lookup_df = pd.DataFrame(dct_lst)
    lookup_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))

if __name__ == "__main__":
    target = sys.argv[1]
    trgFilepath = "/data/mitsuki/out/altorf/denovo/trg/{}/trg.csv".format(target)
    outFilepath = "../data/{}/query_lookup.csv".format(target)
    main(target, trgFilepath, outFilepath)