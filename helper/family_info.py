#!/usr/bin/env python3

import sys
import pandas as pd

def core(target, family):
    baseDirec = "/home/mitsuki/altorf/denovo"
    strainFilepath="{}/data/{}/strain.lst".format(baseDirec, target)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    clusterFilepath="{}/data/{}/cluster.tsv".format(baseDirec, target)
    cluster_df=pd.read_csv(clusterFilepath, sep="\t")
    family_series = cluster_df[cluster_df["family"]==family][strain_lst].dropna(axis=1).iloc[0]
    return family_series

def family_orf(target, family):
    family_series = core(target, family)
    return list(family_series.values)

def family_strain(target, family):
    family_series = core(target, family)
    return list(family_series.index)

if __name__=="__main__":
    target = sys.argv[1]
    family = sys.argv[2]
    print(family_strain(target, family))
    print(family_orf(target, family))
