#!/usr/bin/env python3

import sys
import pandas as pd
from collections import Counter

sys.path.append("../helper")
from phylogeny import NCBIController

def main(target, lcaFilepath, outFilepath):
    nc = NCBIController()
    lca_df = pd.read_csv(lcaFilepath)
    targetGenus_lst = nc.get_descendant(nc.ncbi.get_name_translator([target])[target][0], rank="genus")
    
    lca2genus={}
    for taxid in set(list(lca_df["lca"])):
        if taxid == -1:
            lca2genus[taxid] = -1
        else:
            lineage = nc.get_lineage(taxid)
            if "genus" in lineage.keys() and lineage["genus"] in targetGenus_lst:
                lca2genus[taxid] = lineage["genus"]
            else:
                lca2genus[taxid] = 0
    
    lca_df["genus"] = [lca2genus[lca] for lca in lca_df["lca"]]
    lca_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    baseDirec = "/data/mitsuki/out/altorf/denovo/trg/{}".format(target)
    lcaFilepath = "{}/lca.csv".format(baseDirec)
    outFilepath = "{}/trg.csv".format(baseDirec)
    main(target, lcaFilepath, outFilepath)
