#!/usr/bin/env python3

import sys
import pandas as pd
from ete3 import PhyloTree

sys.path.append("../helper")
from myio import get_strain_lst, get_distance_mat

def main(target, catalogFilepath, nwkFilepath, outFilepath):
    catalogFilepath = "../data/{}/catalog.tsv".format(target)
    catalog_df = pd.read_csv(catalogFilepath, sep="\t")

    nwkFilepath = "../data/{}/cluster.phb".format(target)
    tree = PhyloTree(nwkFilepath)

    #rooting
    outgroup_lst = list(catalog_df[catalog_df["outgroup"] == 1]["genome_id"])
    ancestor = tree.get_common_ancestor(outgroup_lst)
    tree.set_outgroup(ancestor)

    #pruning
    represent_lst = list(catalog_df[catalog_df["represent"] == 1]["genome_id"])
    tree.prune(outgroup_lst + represent_lst, preserve_branch_length=True)

    #output
    tree.write(outfile=outFilepath)
    print("DONE: output {}".format(outFilepath))

if __name__ == "__main__":
    target = sys.argv[1]
    catalogFilepath = "../data/{}/catalog.tsv".format(target)
    nwkFilepath = "../data/{}/cluster.phb".format(target)
    outFilepath = "../data/{}/cluster.format.phb".format(target)
    main(target, catalogFilepath, nwkFilepath, outFilepath)

