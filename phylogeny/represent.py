#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

sys.path.append("../helper")
from myio import *

def main(target, thres, catalogFilepath):
    catalog_df = pd.read_csv(catalogFilepath, sep="\t")

    # filter strain greedy
    distance_mat = get_distance_mat(target, full=True)
    msk = (catalog_df["outgroup"] == 0).values
    for i in range(len(msk)):
        if msk[i]:
            for j in range(i + 1, len(msk)):
                if distance_mat[i][j] < thres:
                    msk[j] = False
    print("DONE: select representative {}/{} strains".format(msk.sum(), len(msk)))

    catalog_df["represent"] = msk.astype(int)
    catalog_df.to_csv(catalogFilepath, index=False, sep="\t")
    print("DONE: update {}".format(catalogFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    thres = float(sys.argv[2])
    catalogFilepath = "../data/{}/catalog.tsv".format(target)
    main(target, thres, catalogFilepath)
