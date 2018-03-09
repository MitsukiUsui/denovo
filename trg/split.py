#!/usr/bin/env python3

import sys
import os
import pandas as pd

sys.path.append("../helper")
from myio import get_cluster_df

def get_header(fp):
    with open(fp, "r") as f:
        header = f.readline().strip()
    return header

target = sys.argv[1]
baseDirec = "/data/mitsuki/out/altorf/denovo/trg/{}".format(target)
resultFilepath = "{}/result.csv".format(baseDirec)
outDirec = "{}/family".format(baseDirec)

header = get_header(resultFilepath)
os.makedirs(outDirec, exist_ok=True)

cluster_df = get_cluster_df(target)
for family in cluster_df["family"]:
    fp = "{}/{}.csv".format(outDirec, family)
    with open(fp, "w") as f:
        f.write("{}\n".format(header))
print("DONE: created {} .csv under {}".format(cluster_df.shape[0], outDirec))

