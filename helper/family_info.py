#!/usr/bin/env python3

import sys
import pandas as pd

sys.path.append("../helper")
from myio import get_strain_lst, get_cluster_df

def core(target, family):
    strain_lst = get_strain_lst(target)
    cluster_df = get_cluster_df(target)
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
