#!/usr/bin/env python3

import sys
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

sys.path.append("../helper")
from myio import get_strain_lst, get_hit_df, get_cluster_df, DistanceManager
from myutil import myinterval

def main(target, strain, outFilepat):

    # summarize query information into list of dictionary
    dct_lst = []
    queryFilepath = "./query/{}/{}.fna".format(target, strain)
    for rec in SeqIO.parse(queryFilepath, "fasta"):
        dct = {
            "family": rec.description.split()[-1],
            "query": rec.id,
            "qlength": len(rec.seq)
        }
        dct_lst.append(dct)

    # extract hit information into myinterval.Interval to prepare for coverage
    hitFilepath = "./result/{}/{}.tsv".format(target, strain)
    hit_df = get_hit_df(hitFilepath)
    query2interval = defaultdict(list)
    for _, row in hit_df.iterrows():
        if row["sstart"] < row["send"]:
            id_ = ";".join([row["sseqid"], str(row["sstart"]-1), str(row["send"]), "0"])
        else:
            id_ = ";".join([row["sseqid"], str(row["send"]-1), str(row["sstart"]), "1"])
        interval = myinterval.Interval(row["qstart"] - 1, row["qend"], id_)
        query2interval[row["qseqid"]].append(interval)

    # calcurate coverage and other information and update list of dict
    dm = DistanceManager(target)
    for dct in dct_lst:
        dct["coverage"] = myinterval.coverage(query2interval[dct["query"]], start=0, end=dct["qlength"])
        dct["distance"] = dm.distance(strain, dct["query"].split('-')[0])
        dct["hit_count"] = len(query2interval[dct["query"]])

        if dct["hit_count"] > 0:
            maxLength = 0
            for interval in query2interval[dct["query"]]:
                if len(interval) > maxLength:
                    maxLength = len(interval)
                    longest = interval.id_
            dct["max_length"] = maxLength
            dct["longest"] = longest

    out_df = pd.DataFrame(dct_lst)
    out_df = out_df[["family", "query", "qlength", "coverage", "distance", "hit_count", "max_length", "longest"]]
    cluster_df = (get_cluster_df(target))[["family", strain]]
    cluster_df["possess"] = cluster_df[strain].notnull().astype(int)
    out_df = pd.merge(out_df, cluster_df[["family", "possess"]], on="family", how="left")
    out_df["traceable"] = ((out_df["qlength"] >= 150) & (out_df["possess"] == 0) &(out_df["coverage"] >= 0.5)).astype(int)
    out_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    outDirec = "./out/{}".format(target)
    os.makedirs(outDirec, exist_ok=True)

    strain_lst = get_strain_lst(target)
    for strain in strain_lst:
        outFilepath = "{}/{}.trace".format(outDirec, strain)
        main(target, strain, outFilepath)
