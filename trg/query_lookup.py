#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO

sys.path.append("../helper")
from myio import get_cluster_df, get_distance_mat, get_strain_lst

def main(target, lcaFilepath, outFilepath):
    full = True

    lca_df = pd.read_csv(lcaFilepath)
    trg_lst = list(lca_df[lca_df["trg"]==1]["family"])
    print("START: create lookup table for {} trg.".format(len(trg_lst), flush=True))

    cluster_df = get_cluster_df(target)
    distance_mat = get_distance_mat(target, full)
    for i in range(distance_mat.shape[0]):
        distance_mat[i][i] = 0
    strain_lst = get_strain_lst(target, full)
    assert len(strain_lst) == distance_mat.shape[0]

    dct_lst = []
    for family in trg_lst:
        dct = {"family": family}

        row = cluster_df[cluster_df["family"]==family].iloc[0]
        mask = row[strain_lst].isnull().values
        for sidx in range(len(strain_lst)): #sbjct index
            x = np.ma.array(distance_mat[sidx], mask=mask)
            qidx=x.argmin() #query index
            orfid = row[strain_lst[qidx]].split(',')[0]
            dct[strain_lst[sidx]] = orfid
        dct_lst.append(dct)

    lookup_df = pd.DataFrame(dct_lst)
    lookup_df = lookup_df[["family"] + strain_lst]
    lookup_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    lcaFilepath = "/data/mitsuki/out/altorf/denovo/trg/{}/lca.csv".format(target)
    outFilepath = "../trg-trace/lookup/{}.csv".format(target)
    main(target, lcaFilepath, outFilepath)
