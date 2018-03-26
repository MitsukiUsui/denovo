#!/usr/bin/env python3

"""
select queryies for mmseqs search against nr
"""

import sys
import os
from Bio import SeqIO

sys.path.append("../helper")
from myio import get_cluster_df, sample_record

def main(target, queryFilepath):
    cluster_df = get_cluster_df(target)
    cluster_df = cluster_df[cluster_df["lineage"] > 1] #dismiss singlton family

    print("START: select {} queries for nr".format(cluster_df.shape[0]), flush=True)
    query_lst = []
    for family in cluster_df["family"]:
        query_lst += sample_record(target, family, "faa", n=1)
    assert len(query_lst) == cluster_df.shape[0]

    os.makedirs(os.path.dirname(os.path.abspath(queryFilepath)), exist_ok=True)
    with open(queryFilepath, 'w') as f:
        for query in query_lst:
            SeqIO.write(query, f, "fasta")
    print("DONE: output {}".format(queryFilepath), flush=True)

if __name__=="__main__":
    target = sys.argv[1]
    queryFilepath = sys.argv[2]
    main(target, queryFilepath)
