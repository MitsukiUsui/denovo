#!/usr/bin/env python3

import sys
import pandas as pd

target=sys.argv[1]
strainFilepath="../data/{}/strain.lst".format(target)
strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

# check result
for strain in strain_lst:
    hitFilepath = "./result/{}/{}.csv".format(target, strain)
    try:
        hit_df=pd.read_csv(hitFilepath)
    except FileNotFoundError:
        print("ERROR: {} does not exist".format(hitFilepath), file = sys.stderr)
        exit(1)
exit(0)
