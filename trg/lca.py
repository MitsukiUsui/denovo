#!/usr/bin/env python3

import sys
import math
import random
import pandas as pd
from ete3 import NCBITaxa

sys.path.append("../helper")
from myio import get_cluster_df
from TaxDbController import TaxDbController
from phylogeny import NCBIController

ncbi = NCBITaxa()
nc = NCBIController()

def get_lca(taxid_lst):
    def get_root(taxid_lst):
        if len(taxid_lst) == 1:  #get_topology() can not handle singleton
            return taxid_lst[0]
        try:
            tree = ncbi.get_topology(taxid_lst)
            root = int(tree.get_tree_root().name)
        except KeyError:
            print("ERROR: get_root() failed", file=sys.stderr)
            return -1
        else:
            return root

    def filterValid(taxid_lst):
        return list(ncbi.get_taxid_translator(taxid_lst).keys())

    taxid_lst = list(set(taxid_lst))
    taxid_lst = filterValid(taxid_lst)

    if len(taxid_lst) == 0:
        print("WARN: no valid taxid found for calculating lca", file=sys.stderr)
        return -1
    else:
        random.shuffle(taxid_lst)
        stp = 100
        for i in range(math.ceil(len(taxid_lst) / stp)):
            end = (i + 1) * stp
            lca = get_root(taxid_lst[:end])
            if lca in [1, 2, 131567]: #root, bacteria, cellular organisms
                    break
        return lca

def main(target, familyDirec, outFilepath):
    cluster_df = get_cluster_df(target)
    family_lst = list(cluster_df["family"])
    print("START: process {} families".format(len(family_lst)), flush=True)

    # assign lca for each family first
    batch = len(family_lst) * 0.01
    border = batch
    dct_lst = []
    for _, family in enumerate(family_lst):
        if _ >= border:
            border += batch
            print(".", end="", flush=True)

        fp = "{}/{}.csv".format(familyDirec, family)
        hit_df = pd.read_csv(fp)
        msk = hit_df["length"] >= (hit_df["qlength"] * 0.2)
        hit_df = hit_df[msk]

        dct = {
            "family" : family,
            "hit_count": hit_df.shape[0],
            "high_count": msk.sum(),
            "query_count": len(set(hit_df["orf_id"])),
            "sbjct_count": len(set(hit_df["accession_version"])),
            "lca": get_lca(hit_df["tax_id"])
        }
        dct_lst.append(dct)
    print()
    lca_df = pd.DataFrame(dct_lst)
    lca_df = lca_df[["family", "hit_count", "high_count", "query_count", "sbjct_count", "lca"]]

    # create lookup table for whether trg or not
    dct_lst = []
    for taxid in set(lca_df["lca"]):
        if taxid == -1:
            dct = {"lca": -1,
                   "trg": False}
        else:
            dct = {"lca": taxid,
                   "trg": nc.is_descendant(child=taxid, parent=nc.encode(target)[0])}
        dct_lst.append(dct)
    lookup_df = pd.DataFrame(dct_lst)
    lookup_df["trg"] = lookup_df["trg"].astype(int)

    # add trg column by merge
    out_df = pd.merge(lca_df, lookup_df, on="lca", how="left")
    out_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))

if __name__ == "__main__":
    target = sys.argv[1]
    baseDirec = "/data/mitsuki/out/altorf/denovo/trg/{}".format(target)
    familyDirec = "{}/family".format(baseDirec)
    outFilepath = "{}/lca.csv".format(baseDirec)
    main(target, familyDirec, outFilepath)
